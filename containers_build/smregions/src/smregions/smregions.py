"""
Contains the command line parsing and the main class of the method
"""

import io
import os
import sys
import logging
from os.path import join, exists
from multiprocessing.pool import Pool

import numpy as np
import pandas as pd
import statsmodels.stats.multitest as mt
from scipy import stats

from smregions import __version__, reference,  load, walker
from smregions.config import file_exists_or_die, file_name
from smregions.executor import ElementExecutor
from smregions.utils import executor_run, loop_logging

logger = logging.getLogger(__name__)


class SMRegions(object):
    """

    Args:
       mutations_file: Mutations input file
       elements_file: Genomic element input file
       regions_of_interest_file: Genemoc elements of interest
       output_folder: Folder where the results will be stored
       config: configuration
       blacklist: File with sample ids (one per line) to remove when loading the input file

    """

    def __init__(self, mutations_file, elements_file, regions_of_interest_file, signature_file, output_file, config):
        logger.debug('Using SMRegions version %s', __version__)

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

        genome_reference_build = self.configuration['reference_genome']
        reference.set_build(genome_reference_build)

        self.cores = self.configuration['cores']
        if self.cores is None:
            self.cores = os.cpu_count()
        logger.debug('Using %s cores', self.cores)

        # Optional parameters
        self.output_file = output_file
        logger.debug('Output: %s', self.output_file)

        self.avoid_parallel = False

        s = io.BytesIO()
        self.configuration.write(s)
        logger.debug('Configuration used:\n' + s.getvalue().decode())
        s.close()

        np.random.seed(self.configuration['seed'])

    def run(self):
        """
        Run the analysis.
        """

        # Load mutations mapping
        mutations, elements, regions_of_interest = load.load_and_map_variants(self.mutations_file,
                                                                        self.elements_file,
                                                                        self.regions_of_interest_file)

        # Load signatures
        if self.signature_file is None:
            signature = None
        else:
            logger.debug('Loading signature')
            signature = pd.read_pickle(self.signature_file)['probabilities']

        # Create one executor per element
        element_executors = [ElementExecutor(element_id, muts, elements[element_id],
                                             regions_of_interest[element_id],
                                             signature, self.configuration, np.random.randint(0, 2**32-1))
                             for element_id, muts in sorted(mutations.items())
                             if len(muts) >= self.configuration['muts_min']]

        # Sort executors to compute first the ones that have more mutations
        element_executors = sorted(element_executors, key=lambda e: -len(e.muts))

        # Run the executors
        with Pool(self.cores) as pool:
            results = {}
            logger.info("Computing SMRegions")
            map_func = map if self.avoid_parallel else pool.imap
            for executor in loop_logging(map_func(executor_run, element_executors), size=len(element_executors), step=6*self.cores):
                results[executor.name] = executor.result

            # Flatten partitions
            partitions = list(walker.flatten_partitions(results))

            i = 0
            while len(partitions) > 0 or i == 0:

                i += 1
                logger.info("Parallel sampling. Iteration %d, genes %d, partitions %d", i, len(set([n for n,p,r,s in partitions])), len(partitions))

                # Pending sampling execution
                for _ in loop_logging(map_func(walker.compute_sampling, partitions), size=len(partitions), step=1):
                    continue

        # Compute p-values
        logger.info("Computing p-values")
        list_results = []

        for result in results.values():
            total_mut = result['nmuts']
            for region_name, region_counts in result['in_reg_counts'].items():
                observed, simulated = region_counts
                mean_simulated = simulated / self.configuration['sampling']
                a = observed
                b = total_mut - observed
                c = mean_simulated
                d = total_mut - mean_simulated
                if a > 0:
                    u,p_value = stats.power_divergence(f_obs=[a, b], f_exp=[c, d], lambda_="log-likelihood")
                    list_results.append([region_name.split(";")[0], region_name.split(";")[1], total_mut,
                                         observed, mean_simulated, u, p_value])

        if results == {}:
            logger.warning("Empty resutls, possible reason: no mutation from the dataset can be mapped to the provided regions.")
            sys.exit(0)

        # Sort and store results
        logger.info("Storing results")
        df_results = pd.DataFrame(list_results,
                                  columns=["REGION", "HUGO_SYMBOL", "TOTAL_MUTS_GENE", "OBSERVED_REGION", "MEAN_SIMULATED",
                                           "U", "P_VALUE"])
        
        # Correct the p-value
        if len(df_results) > 0:
            p = df_results["P_VALUE"].values
            mask = np.isfinite(p)
            pval_corrected = np.full(p.shape, np.nan)
            pval_corrected[mask] = mt.multipletests(p[mask], method='fdr_bh')[1]
            df_results["Q_VALUE"] = pval_corrected
        else:
            df_results["Q_VALUE"] = df_results["P_VALUE"]

        df_results.to_csv(self.output_file, sep="\t", index=False, compression="gzip")
        logger.info("Done")
