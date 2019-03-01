import gzip
import os
import subprocess
import sys
import tempfile

import click
import pandas as pd


# global variables

folder = '/deconstructsig'

DRIVERS_PATH = os.path.join(folder, 'output_pass_drivers_01.csv')
tmpdir = tempfile.mkdtemp()
MUTS_PATH = os.path.join(tmpdir, 'mutations_snv.csv')


# run functions

@click.command()
@click.option('--input_file', type=click.Path(), help='MAF file')
@click.option('--weights', type=click.Path(), help='output table with signature weights')
@click.option('--build', type=str, help='reference genome build')
def run_deconstruct(input_file, weights, build):

    """
    Run deconstructSigs.r against COSMIC signatures.
    Returns weights of active COSMIC signatures for each sample.
    """

    weights_r = os.path.splitext(weights)[0]
    cohort = os.path.basename(input_file).split('.')[0]

    # read input
    muts = pd.read_csv(input_file, sep='\t', compression="gzip")

    # filter out non SNVs
    muts = muts[(muts['REF'] != '-') & (muts['REF'].isin(list('ACGT')))]
    muts = muts[(muts['ALT'] != '-') & (muts['ALT'].isin(list('ACGT')))]

    # filter out mutations in driver genes
    print(cohort)
    df_filtered = pd.read_csv(DRIVERS_PATH, sep="\t")
    df_filtered = df_filtered[df_filtered["COHORT"] == cohort]

    # if df_filtered.shape[0] > 0:
    #     list_drivers = df_filtered["SYMBOL"].unique()
    #     muts = muts[~(muts["GENE"].isin(list_drivers))]

    # keep mutations in a temporal file
    muts.to_csv(MUTS_PATH, sep='\t', index=False)

    # run deconstructSigs
    command = 'Rscript ' + os.path.join(os.path.abspath(os.path.dirname(__file__)), 'deconstructSigs.r') + \
              ' ' + MUTS_PATH + ' ' + build + ' ' + weights_r

    print(command)
    r = subprocess.run(command, shell=True, executable='/bin/bash')
    if r.returncode != 0:
        print("ERROR: deconstructSigs.r failed")
        sys.exit(r.returncode)

    # exit without error if no signatures_weight.csv exist
    if not os.path.exists(weights_r):
        print("WARNING! no sample got enough variants to undergo fitting")
        sys.exit(0)

    # Compress the file
    with open(weights_r, 'r') as fi, gzip.open(weights, 'wt') as fo:
        fo.write(fi.read())


if __name__ == '__main__':

    run_deconstruct()
