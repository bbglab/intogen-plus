"""
Contains the command line parsing and the main class of the method
"""

import io
import os
import sys
import csv
import logging
import numpy as np
from scipy import stats

from os.path import join, exists
from collections import defaultdict
from multiprocessing.pool import Pool
from logging.config import dictConfig

import click
import pandas as pd
import bgdata
from bgpack import packer, reader as bgpack_reader

from simReg import __version__, __logger_name__
from simReg.config import load_configuration, file_exists_or_die, file_name
from simReg.executors.element import ElementExecutor
from simReg.load import load_and_map_variants
from simReg.signature import change_ref_build
from simReg.store import store_tsv, store_png, store_html
from simReg.utils import executor_run, loop_logging
from simReg.walker import flatten_partitions, compute_sampling

logger = logging.getLogger(__logger_name__)

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


class OncodriveFML(object):
    """

    Args:
       mutations_file: Mutations input file (see :mod:`~simReg.load` for details)
       elements_file: Genomic element input file (see :mod:`~simReg.load` for details)
       output_folder: Folder where the results will be stored
       config: configuration (see :ref:`configuration <project configuration>`)
       blacklist: File with sample ids (one per line) to remove when loading the input file
       generate_pickle (bool): flag to only generate pickle files and exit

    """

    def __init__(self, mutations_file, elements_file, regions_of_interest_file, signature_file, output_folder, config, blacklist, generate_pickle):
        logger.debug('Using OncodriveFML version %s', __version__)

        # Required parameters
        self.mutations_file = file_exists_or_die(mutations_file)
        logger.debug('Mutations file: %s', self.mutations_file)
        self.elements_file = file_exists_or_die(elements_file)
        logger.debug('Elements file: %s', self.elements_file)
        self.regions_of_interest_file = file_exists_or_die(regions_of_interest_file)
        logger.debug('Elements of interest file: %s', self.regions_of_interest_file)
        self.signature_file = signature_file
        logger.debug('Signature file: %s', self.signature_file)
        self.configuration = config
        self.blacklist = blacklist
        if self.blacklist is not None:
            logger.debug('Blacklist file: %s', self.blacklist)
            logger.debug('Using a blacklist causes some pickle files not to be saved/loaded')
        self.generate_pickle = generate_pickle

        genome_reference_build = self.configuration['genome']['build']
        change_ref_build(genome_reference_build)

        self.cores = self.configuration['settings']['cores']
        if self.cores is None:
            self.cores = os.cpu_count()
        logger.debug('Using %s cores', self.cores)

        # Optional parameters
        self.output_folder = output_folder
        self.output_file_prefix = join(self.output_folder, file_name(self.mutations_file) + '-simReg')
        logger.debug('Output: %s', self.output_folder)

        # Output parameters
        self.mutations = None
        self.elements = None
        self.regions_of_interest = None
        self.signature = None

        self.avoid_parallel = False

        s = io.BytesIO()
        self.configuration.write(s)
        logger.debug('Configuration used:\n' + s.getvalue().decode())
        s.close()

    def run(self):
        """
        Run the OncodriveFML analysis.
        """

        # Load mutations mapping
        mutations_data, self.elements, self.regions_of_interest = load_and_map_variants(self.mutations_file,
                                                              self.elements_file,
                                                              self.regions_of_interest_file,
                                                              blacklist=self.blacklist,
                                                              save_pickle=self.generate_pickle)

        self.mutations = mutations_data['data']

        if self.generate_pickle:
            logger.info('Pickles generated. Exiting')
            return

        # Load signatures
        if self.signature_file is None:
            self.signature = None
        else:
            logging.debug('Loading signature')
            self.signature = pd.read_pickle(self.signature_file)['probabilities']

        # self.vep_reader = bgpack_reader.load(bgdata.get_path('bgvep', 'GRCh37', 'most_severe'))

        # Create one executor per element
        element_executors = [ElementExecutor(element_id, muts, self.elements[element_id],
                                             self.regions_of_interest[element_id],
                                             None,
                                             self.signature, self.configuration) for
                             element_id, muts in self.mutations.items()]

        # Remove elements that don't have mutations to compute
        element_executors = [e for e in element_executors if len(e.muts) > 0]

        # Sort executors to compute first the ones that have more mutations
        element_executors = sorted(element_executors, key=lambda e: -len(e.muts))

        # Run the executors
        with Pool(self.cores) as pool:
            results = {}
            logger.info("Computing OncodriveFML")
            map_func = map if self.avoid_parallel else pool.imap
            for executor in loop_logging(map_func(executor_run, element_executors), size=len(element_executors), step=6*self.cores):
                if executor.result['muts_count'] > 0:
                    results[executor.name] = executor.result

            # Flatten partitions
            partitions = list(flatten_partitions(results))

            i = 0
            while len(partitions) > 0 or i == 0:

                i += 1
                logger.info("Parallel sampling. Iteration %d, genes %d, partitions %d", i, len(set([n for n,p,r in partitions])), len(partitions))

                # Pending sampling execution
                for _ in loop_logging(map_func(compute_sampling, partitions), size=len(partitions), step=1):
                    continue

            # Compute p-values
            logger.info("Compute p-values")
            list_results = []

            for result in results.values():

                total_mut = result['muts_count']

                for region_name, region_counts in result['in_reg_counts'].items():
                    observed = region_counts.pop(0)
                    mean_simulated = np.mean(region_counts)
                    a = observed
                    b = total_mut - observed
                    c = mean_simulated
                    d = total_mut - mean_simulated
                    if a > 0:
                        u,p_value = stats.power_divergence(f_obs=[a, b], f_exp=[c, d], lambda_="log-likelihood")
                        list_results.append([region_name.split(";")[0],region_name.split(";")[1],total_mut,observed,mean_simulated,u,p_value])




        if results == {}:
            logger.warning("Empty resutls, possible reason: no mutation from the dataset can be mapped to the provided regions.")
            sys.exit(0)

        # Sort and store results
        logger.info("Storing results")
        if not exists(self.output_folder):
            os.makedirs(self.output_folder, exist_ok=True)
        df_results = pd.DataFrame(list_results,
                                  columns=["DEGRON","SYMBOL", "TOTAL_MUTS_GENE", "OBSERVED_REGION", "MEAN_SIMULATED", "U",
                                           "P_VALUE"])
        results_file = self.output_file_prefix + ".tsv"
        df_results.to_csv(results_file,sep="\t",index=False)

        logger.info("Creating figures")
        store_png(results_file, self.output_file_prefix + ".png")
        #store_html(result_file, self.output_file_prefix + ".html")

        logger.info("Done")


def main(mutations_file, elements_file, regions_file, signature_file, output_folder, config_file, samples_blacklist, generate_pickle, config_override_dict=None):

    configuration = load_configuration(config_file, override=config_override_dict)

    output_folder = file_name(elements_file) if output_folder is None else output_folder
    output_file = join(output_folder, file_name(mutations_file) + '-simReg.tsv')
    # Skip if done
    if exists(output_file):
        store_png(output_file, output_file + ".png")
        logging.warning("Already calculated at '{}'".format(output_file))
        return
    else:
        if not exists(output_folder):
            os.makedirs(output_folder, exist_ok=True)
        if 'logging' in configuration:
            if configuration['logging'].get('handlers', {}).get('file', None) is not None:
                logging_file = join(output_folder, file_name(mutations_file) + '__log.txt')
                configuration['logging']['handlers']['file']['filename'] = logging_file
            dictConfig(configuration['logging'])

    analysis = OncodriveFML(mutations_file, elements_file, regions_file, signature_file, output_folder, configuration,
                            samples_blacklist, generate_pickle)

    logger.info('Running analysis')
    # Run the analysis
    analysis.run()


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--input', 'mutations_file', type=click.Path(exists=True), help='Variants file', metavar='MUTATIONS_FILE',required=True)
@click.option('-e', '--elements', 'elements_file', type=click.Path(exists=True), metavar='ELEMENTS_FILE', help='Genomic elements to analyse', required=True)
@click.option('-r', '--regions', 'regions_file', type=click.Path(exists=True), metavar='SEGMENTS_FILE', help='Genomic segments of interest', required=True)
@click.option('-sg', '--signature', 'signature_file', type=click.Path(exists=True), metavar='SIGNATURE_FILE', help='Signature file', default=None)
@click.option('-o', '--output', 'output_folder', type=click.Path(), metavar='OUTPUT_FOLDER', help="Output folder. Default to regions file name without extensions.", default=None)
@click.option('-c', '--configuration', 'config_file', default=None, type=click.Path(exists=True), metavar='CONFIG_FILE', help="Configuration file. Default to 'oncodrivefml_v2.conf' in the current folder if exists or to ~/.bbglab/oncodrivefml_v2.conf if not.")
@click.option('--samples-blacklist', default=None, type=click.Path(exists=True), metavar='SAMPLES_BLACKLIST', help="Remove these samples when loading the input file.")
@click.option('--generate-pickle', help="Run OncodriveFML to generate pickle files that could speed up future executions and exit.", is_flag=True)
@click.option('--debug', help="Show more progress details", is_flag=True)
@click.version_option(version=__version__)
def cmdline(mutations_file, elements_file, regions_file, signature_file, output_folder, config_file, samples_blacklist, generate_pickle, debug):
    """
    Run OncodriveFML on the genomic regions in ELEMENTS FILE
    using the mutations in MUTATIONS FILE.

    """
    dd = lambda: defaultdict(dd)
    override_config = dd()

    if debug:
        override_config['logging']['handlers']['console']['level'] = 'DEBUG'
    else:
        override_config['logging']['handlers']['console']['level'] = 'INFO'

    main(mutations_file, elements_file, regions_file, signature_file, output_folder, config_file, samples_blacklist, generate_pickle, override_config)


if __name__ == "__main__":
    cmdline()
