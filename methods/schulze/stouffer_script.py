import os

import click as click
import pandas as pd
import numpy as np
from scipy.stats import combine_pvalues, uniform
import statsmodels.sandbox.stats.multicomp as multicomp


DEFAULT_METHODS = ['oncodrivefml', 'oncodriveclust', 'oncodriveomega', 'hotmapssignature', 'mutsigcv']


def parse_optimized_weights(path_weights):
    cap = lambda a: a[:-2]
    df = pd.read_csv(path_weights, sep='\t', compression="gzip")
    del df['Objective_Function']
    dict_weight = df.to_dict()
    return {cap(k): v[0] for k, v in dict_weight.items()}


def retrieve_ranking(df, path_ranking):
    """
    df: dataframe with p-values per method
    tumor_type: str
    """

    ranking_dict = pd.read_csv(path_ranking,
                       sep='\t',
                       usecols=['SYMBOL', 'RANKING', 'Median_Ranking', 'Total_Bidders', 'All_Bidders'],
                       low_memory=False,
                               compression="gzip"
                       )
    cols = ['SYMBOL']
    return pd.merge(left=df, right=ranking_dict, left_on=cols, right_on=cols, how="left")


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


def trim_nans(pvals):
    """
    provides reduced list of pvalues removing nans
    """

    nan_mask = np.isnan(pvals)
    reduced_pvals = pvals[~nan_mask]
    return reduced_pvals, nan_mask


def truncate(pvals, threshold=1e-16):

    mask = (pvals < threshold)
    pvals[mask] = threshold
    return pvals


def trimmed_stouffer_w(pvals, weights):
    """
    conducts stouffer_w where pvals and weights are clean from nans
    """

    reduced_pvals, nan_mask = trim_nans(pvals)
    reduced_weights = weights[~nan_mask]
    return stouffer_w(truncate(reduced_pvals), weights=reduced_weights)
def load_cgc():
    '''
    Loads the CGC set and returns a set of CGC genes
    :return: set of cgc gene names
    '''
    df_cgc = pd.read_csv("//workspace/projects/intogen_2017/data/latest/CGC_set.tsv",sep="\t")
    return set(df_cgc["Gene Symbol"].values)
def set_qvalue_cgc(row,qvalues_cgc,cgc_set):
    '''
    Set the CGC qvalue to rows that are CGC genes
    :param row:
    :return: the corrected qvalue or nan
    '''
    i = 0
    if row["SYMBOL"].isin(cgc_set):
        i
def combine_pvals(df, path_weights):

    weight_dict = parse_optimized_weights(path_weights)
    weights = np.array([weight_dict[m] for m in DEFAULT_METHODS])
    func = lambda x: trimmed_stouffer_w(x, weights)
    df['PVALUE_' + 'stouffer_w'] = df[['PVALUE_' + m for m in DEFAULT_METHODS]].apply(func, axis=1)
    df['QVALUE_' + 'stouffer_w'] = multicomp.multipletests(df['PVALUE_' + 'stouffer_w'].values, method='fdr_bh')[1]
    cgc_set = load_cgc()
    df_cgc = df[df["SYMBOL"].isin(cgc_set)].copy()
    pvalues_cgc = df_cgc['PVALUE_' + 'stouffer_w'].values
    qvalues_cgc = multicomp.multipletests(pvalues_cgc, method='fdr_bh')[1]
    df_cgc['QVALUE_CGC_' + 'stouffer_w'] = qvalues_cgc
    df_final = pd.merge(left=df,right=df_cgc[["SYMBOL","QVALUE_CGC_stouffer_w"]],left_on="SYMBOL",right_on=["SYMBOL"],how="left")
    return df_final


def partial_correction(df, fml_data):

    dh = pd.merge(left=df, right=fml_data[['SYMBOL', 'Q_VALUE']], left_on=['SYMBOL'], right_on=['SYMBOL'], how="left")
    c = dh['Q_VALUE'].values
    mask = ~np.isnan(c)
    a = dh['PVALUE_' + 'stouffer_w'].values

    if len(a[mask]) == 0:
        raise RuntimeError("No data after filtering NaN OncodriveFML q-values")

    c[mask] = multicomp.multipletests(a[mask], method='fdr_bh')[1]
    dh['QVALUE_' + 'stouffer_w'] = c
    del dh['Q_VALUE']
    return dh


def combine_from_tumor(df, path_to_output, path_fml):
    fml_data = pd.read_csv(path_fml, sep='\t', compression="gzip")
    dh = partial_correction(df, fml_data)
    dh.to_csv(path_to_output, sep='\t', index=False, compression="gzip")


@click.command()
@click.option('--input_path', type=click.Path(exists=True), help="Path to input dataframe with SYMBOL and p-value for each method", required=True)
@click.option('--output_path', type=click.Path(), help="Path to output dataframe", required=True)
@click.option('--path_rankings', type=click.Path(), help="Path to dataframe produced by the voting system", required=True)
@click.option('--path_weights', type=click.Path(), help="Path to dataframe with weights", required=True)
@click.option('--path_fml', type=click.Path(), help="Path to OncodriveFML results folder", required=True)
def run_stouffer_script(input_path, output_path, path_rankings, path_weights, path_fml):

    df = pd.read_csv(input_path, sep='\t', compression="gzip")
    dg = retrieve_ranking(df, path_rankings)
    dh = combine_pvals(dg, path_weights)
    combine_from_tumor(dh, output_path, path_fml)


if __name__ == '__main__':
    run_stouffer_script()

