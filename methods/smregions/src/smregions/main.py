"""
Contains the command line parsing
"""

import os
from os import path

import click
import bglogs

from smregions import __version__, config
from smregions.smdeg import SMDeg


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def main(mutations_file, elements_file, regions_file, signature_file, output_folder, config_file):

    output_folder = config.file_name(elements_file) if output_folder is None else output_folder
    output_file = path.join(output_folder, config.file_name(mutations_file) + '-smregions.tsv')
    # Skip if done
    if path.exists(output_file):
        bglogs.warning("Already calculated at '{}'".format(output_file))
        return
    else:
        if not path.exists(output_folder):
            os.makedirs(output_folder, exist_ok=True)

    configuration = config.load(config_file)

    analysis = SMDeg(mutations_file, elements_file, regions_file, signature_file, output_folder, configuration)

    bglogs.info('Running analysis')
    # Run the analysis
    analysis.run()


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-m', '--muts', 'mutations_file', type=click.Path(exists=True), help='Variants file', metavar='MUTATIONS_FILE',required=True)
@click.option('-e', '--elements', 'elements_file', type=click.Path(exists=True), metavar='ELEMENTS_FILE', help='Genomic elements to analyse', required=True)
@click.option('-r', '--regions', 'regions_file', type=click.Path(exists=True), metavar='REGIONS_FILE', help='Genomic regions of interest', required=True)
@click.option('-s', '--signature', 'signature_file', type=click.Path(exists=True), metavar='SIGNATURE_FILE', help='Signature file. Default equial probabilities', default=None)
@click.option('-o', '--output', 'output_folder', type=click.Path(), metavar='OUTPUT_FOLDER', help="Output folder. Default to regions file name without extensions.", default=None)
@click.option('-c', '--configuration', 'config_file', default=None, type=click.Path(exists=True), metavar='CONFIG_FILE', help="Configuration file. Default to 'smregions.conf' in the current folder if exists or to ~/.bbglab/simreg.conf if not.")
@click.option('--debug', help="Show more progress details", is_flag=True)
@click.version_option(version=__version__)
def cmdline(mutations_file, elements_file, regions_file, signature_file, output_folder, config_file, debug):
    """
    Run SMDeg on the genomic regions in ELEMENTS FILE and the regions of interest REGIONS_FILE
    using the mutations in MUTATIONS FILE.

    """
    bglogs.configure(debug=True if debug else False)

    main(mutations_file, elements_file, regions_file, signature_file, output_folder, config_file)


if __name__ == "__main__":
    cmdline()
