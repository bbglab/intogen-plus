import os

import click as click
import pandas as pd
import numpy as np
from scipy.stats import combine_pvalues, uniform
import statsmodels.sandbox.stats.multicomp as multicomp

# INTOGEN_PATH = '/workspace/projects/intogen/intogen4/'
# path_ranking = os.path.join(INTOGEN_PATH, 'scripts', 'data', 'results', 'optimization', 'ranking')
# path_weights = os.path.join(INTOGEN_PATH, 'scripts', 'data', 'results', 'optimization', 'weights')
# path_fml = os.path.join(INTOGEN_PATH, 'runs', 'intogen4_20170614', 'oncodrivefml')

tumor_types = ['UCS', 'LUAD', 'LIHC', 'THYM', 'DLBC', 'PAAD', 'KICH', 'PCPG', 'SKCM', 'BRCA',
               'UCEC', 'ACC', 'COADREAD', 'CESC', 'TGCT', 'LAML', 'MESO', 'COAD', 'CHOL', 'KIRC',
               'HNSC', 'LUSC', 'UVM', 'ESCA', 'BLCA', 'SARC', 'READ', 'OV', 'PRAD', 'THCA', 'LGG',
               'KIRP', 'GBM', 'STAD']

default_methods_list = ['oncodrivefml', 'oncodriveclust', 'oncodriveomega', 'hotmapssignature', 'mutsigcv']


def parse_optimized_weights(path_weights):

    cap = lambda a: a[:-2]
    df = pd.read_csv(path_weights, sep='\t')
    del df['Objective_Function']
    dict_weight = df.to_dict()
    return {cap(k): v[0] for k, v in dict_weight.items()}

def create_ranking_dict(tumor_type, path_ranking):

    """
    Creates a dictionary of dataframe labelled by tumor type
    """

    ranking_dict = {}
    try:
        ranking_dict[tumor_type] = pd.read_csv(path_ranking, sep='\t',
                                               usecols=['Cancer_Type', 'SYMBOL', 'RANKING', 'Median_Ranking', 'Total_Bidders', 'All_Bidders'],
                                               low_memory=False)
    except:
        print('{0} was not successful 1'.format(tumor_type))
    return ranking_dict


def retrieve_ranking(df, tumor_type, path_ranking):

    """
    df: dataframe with p-values per method
    tumor_type: str
    """

    ranking_dict = create_ranking_dict(tumor_type, path_ranking)
    try:
        dh = ranking_dict[tumor_type]
        cols = ['SYMBOL']
        return pd.merge(left=df, right=dh, left_on=cols, right_on=cols, how="left")
    except:
        print('{0} was not successful 2'.format(tumor_type))


def stouffer_w(pvals, weights=None):

    return combine_pvalues(pvals, method='stouffer', weights=weights)[1]


def impute(pvals):
    """
    impute array-like instance with uniform [0,1] distribution
    """

    mask1 = np.isnan(pvals)
    mask2 = (pvals == 1)
    mask = mask1 | mask2
    pvals[mask] = uniform.rvs(size=len(pvals[mask]))
    return pvals


def combine_pvals(df, path_weights):

    weight_dict = parse_optimized_weights(path_weights)
    weights = np.array([weight_dict[m] for m in default_methods_list])
    func = lambda x: stouffer_w(impute(x), weights=weights)
    df['PVALUE_' + 'stouffer_w'] = df[['PVALUE_' + m for m in default_methods_list]].apply(func, axis=1)
    df['QVALUE_' + 'stouffer_w'] = multicomp.multipletests(df['PVALUE_' + 'stouffer_w'].values, method='fdr_bh')[1]
    return df


def partial_correction(df, fml_data):

    dh = pd.merge(left=df, right=fml_data[['SYMBOL', 'Q_VALUE']], left_on=['SYMBOL'], right_on=['SYMBOL'], how="left")
    c = dh['Q_VALUE'].values
    mask = ~np.isnan(c)
    a = dh['PVALUE_' + 'stouffer_w'].values
    c[mask] = multicomp.multipletests(a[mask], method='fdr_bh')[1]
    dh['QVALUE_' + 'stouffer_w'] = c
    del dh['Q_VALUE']
    return dh


def combine_from_tumor(df, path_to_output, tumor_type, path_fml):

    path_fml = os.path.join(path_fml, '{0}.out.gz'.format(tumor_type))
    fml_data = pd.read_csv(path_fml, sep='\t')
    dh = partial_correction(df, fml_data)
    dh.to_csv(path_to_output, sep='\t', index=False)


@click.command()
@click.option('--input_path', type=click.Path(exists=True),
              help="Path to input dataframe with SYMBOL and p-value for each method", required=True)
@click.option('--output_path', type=click.Path(),
              help="Path to output dataframe", required=True)
@click.option('--tumor_type', type=click.Path(),
              help="Tumor type ID", required=True)
@click.option('--path_ranking', type=click.Path(),
              help="Path to dataframe produced by the voting system", required=True)
@click.option('--path_weights', type=click.Path(),
              help="Path to dataframe with weights", required=True)
@click.option('--path_fml', type=click.Path(),
              help="Path to OncodriveFML results folder", required=True)

def run_stouffer_script(input_path, output_path, tumor_type, path_ranking, path_weights, path_fml):

    df = pd.read_csv(input_path, sep='\t')
    dg = retrieve_ranking(df, tumor_type, path_ranking)
    dh = combine_pvals(dg, path_weights)
    combine_from_tumor(dh, output_path, tumor_type, path_fml)

if __name__ == '__main__':

    run_stouffer_script()

