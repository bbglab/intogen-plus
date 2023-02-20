import os
import gzip
import pickle
import json
import click
from functools import partial
from tqdm import tqdm

import numpy as np
import pandas as pd

from utils import complementary


# load 96-array with triplet counts
trinuc_subs =  [
    'ACA>A', 'ACA>G', 'ACA>T', 'ACC>A', 'ACC>G', 'ACC>T', 'ACG>A', 'ACG>G', 'ACG>T', 'ACT>A',
    'ACT>G', 'ACT>T', 'ATA>A', 'ATA>C', 'ATA>G', 'ATC>A', 'ATC>C', 'ATC>G', 'ATG>A', 'ATG>C',
    'ATG>G', 'ATT>A', 'ATT>C', 'ATT>G', 'CCA>A', 'CCA>G', 'CCA>T', 'CCC>A', 'CCC>G', 'CCC>T', 
    'CCG>A', 'CCG>G', 'CCG>T', 'CCT>A', 'CCT>G', 'CCT>T', 'CTA>A', 'CTA>C', 'CTA>G', 'CTC>A', 
    'CTC>C', 'CTC>G', 'CTG>A', 'CTG>C', 'CTG>G', 'CTT>A', 'CTT>C', 'CTT>G', 'GCA>A', 'GCA>G', 
    'GCA>T', 'GCC>A', 'GCC>G', 'GCC>T', 'GCG>A', 'GCG>G', 'GCG>T', 'GCT>A', 'GCT>G', 'GCT>T', 
    'GTA>A', 'GTA>C', 'GTA>G', 'GTC>A', 'GTC>C', 'GTC>G', 'GTG>A', 'GTG>C', 'GTG>G', 'GTT>A', 
    'GTT>C', 'GTT>G', 'TCA>A', 'TCA>G', 'TCA>T', 'TCC>A', 'TCC>G', 'TCC>T', 'TCG>A', 'TCG>G', 
    'TCG>T', 'TCT>A', 'TCT>G', 'TCT>T', 'TTA>A', 'TTA>C', 'TTA>G', 'TTC>A', 'TTC>C', 'TTC>G', 
    'TTG>A', 'TTG>C', 'TTG>G', 'TTT>A', 'TTT>C', 'TTT>G']

with open('/mutrate/cds_triplet_abundance.json', 'rt') as f:
    tricount_exome = json.load(f)


with open('/mutrate/genome_triplet_abundance.json', 'rt') as f:
    tricount_genome = json.load(f)


class dNdSOut:

    def __init__(self, annotmuts_fn):

        self.annotmuts = pd.read_csv(annotmuts_fn, sep='\t')

    @staticmethod
    def pyr_context(row):
        """pyr-context associated to row in maf table"""

        context = row['ref3_cod'] + '>' + row['mut_cod']
        if context[1] not in list('CT'):
            context = complementary(context)
        return context

    def catalogue(self):
        """96 channel count matrix"""

        df = self.annotmuts[(self.annotmuts['ref'].isin(list('ACGT'))) & (self.annotmuts['mut'].isin(list(
            'ACGT')))].copy()
        df.loc[:, 'context'] = df.apply(self.pyr_context, axis=1)
        catalogue = pd.crosstab(df.sampleID, df.context)

        # check if we are missing some context, if so add it with a 0-filled column.
        if catalogue.shape[1] < 96:
            missing_context = set(list(trinuc_subs)).difference(set(catalogue.columns))
            for context in missing_context:
                catalogue.loc[:, context] = 0
            
            catalogue = catalogue.reindex(sorted(catalogue.columns), axis=1)

        return catalogue


def cohort_profile(annotmuts_fn):
    """
    :param gene: gene symbol
    :param annotmuts: dndscv annotmuts filename
    :return: dictionary of mutrates per gene
    """

    dndscv = dNdSOut(annotmuts_fn)
    catalogue = dndscv.catalogue()
    res = np.einsum('ij->j', catalogue.values)
    print(res)
    return res


def cohort_norm_profile(profile, scope):
    """
    :param profile: 96 channel array of counts
    :param scope: 'exome' or 'genome'
    :return: 96 channel array of frequencies
    """

    if scope == 'exome':
        res = np.array(profile) / np.array(tricount_exome)
    elif scope == 'genome':
        res = np.array(profile) / np.array(tricount_genome)
    return res / np.sum(res)


@click.group()
def cli():
    pass


@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--annotmuts', 'annotmuts_fn', type=click.Path(), help='path to dndsout$annotmuts table')
@click.option('--scope', type=str, help='either "genome" or "exome"')
@click.option('--outfolder', required=True, type=click.Path(), help='path to output')
def mutrate(annotmuts_fn, scope, outfolder):

    profile = cohort_profile(annotmuts_fn)
    norm_profile = list(cohort_norm_profile(profile, scope))

    label = os.path.basename(annotmuts_fn).split('.')[0]
    outfile = os.path.join(outfolder, f'{label}.mutrate.json')
    with open(outfile, 'wt') as f:
        json.dump(norm_profile, f)


if __name__ == '__main__':

    """
    Example:

    python compute_mutrate.py mutrate --annotmuts /workspace/datasets/intogen/oriolRun/ \
    clonal_hematopoiesis_run_20201015_DBS_out/intogen_20200929/debug/dndscv/ \
    OTHER_WXS_TCGA_FULL.dndscv_annotmuts.tsv.gz --scope exome --outfolder ./

    """

    cli()
