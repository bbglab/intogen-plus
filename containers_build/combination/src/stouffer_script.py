import os

import click as click
import pandas as pd
import numpy as np
from functools import reduce
from scipy.stats import combine_pvalues
from statsmodels.stats.multitest import multipletests

from combination import custom_combination

from configobj import ConfigObj


# Global variables

CODE_DIR = os.path.join(
    os.path.dirname(
        os.path.abspath(__file__)
    )
)

configs_file = os.path.join(CODE_DIR, 'config', 'intogen_qc.cfg')
configs = ConfigObj(configs_file)
METHODS = configs['methods'].keys()


def parse_optimized_weights(path_weights):
    '''
    Parse the dataframe with the weights
    :param path_weights:
    :return: a dictionary with the weights of each method
    '''
    df = pd.read_csv(path_weights, sep='\t', compression="gzip")
    del df['Objective_Function']
    dict_weight = df.to_dict()
    return {k: v[0] for k, v in dict_weight.items()}


def retrieve_ranking(df, path_ranking):
    '''
    Retrieve the ranking of the genes
    :param df: dataframe of input
    :param path_ranking: path of the ranked genes
    :return:
    '''

    df_ranked = pd.read_csv(path_ranking,
                       sep='\t',
                       low_memory=False,
                               compression="gzip"
                       )
    cols = ['SYMBOL']
    return pd.merge(left=df, right=df_ranked, left_on=cols, right_on=cols, how="left")


def stouffer_w(pvals, weights=None):
    '''
    Compute the weighted stouffer p-value
    :param pvals:
    :param weights:
    :return:
    '''

    return combine_pvalues(pvals, method='stouffer', weights=weights)[1]


def trim_nans(pvals):
    '''
    exclude nans for the trimming
    :param pvals:
    :return:
    '''
    nan_mask = np.isnan(pvals) # used to exclude nan p-values from combination
    anti_one_mask = pvals[pvals < 1.] # used to exclude p-values 1.0 from combination
    mask = ~nan_mask & anti_one_mask
    return pvals[mask], mask


def truncate(pvals, threshold=1e-16):
    '''
    Truncate p-values to a minimum value of threshold
    :param pvals:
    :param threshold:
    :return:
    '''
    mask = (pvals < threshold)
    pvals[mask] = threshold
    return pvals


def trimmed_stouffer_w(pvals, weights):
    '''
    Conducts stouffer_w where pvals and weights are clean from nans
    :param pvals:
    :param weights:
    :return:
    '''
    reduced_pvals, mask = trim_nans(pvals)
    reduced_weights = weights[mask]
    return stouffer_w(truncate(reduced_pvals), weights=reduced_weights)


def load_cgc():
    '''
    Loads the CGC set and returns a set of CGC genes
    :return: set of cgc gene names
    '''
    df_cgc = pd.read_csv(os.path.join(os.environ['INTOGEN_DATASETS'], "combination",'cgc', "cancer_gene_census_parsed.tsv"), sep="\t")
    return set(df_cgc["Gene Symbol"].values)


def combine_pvals(df, path_weights, methods):
    '''
    Create the combined p-values of stouffer for all genes and another independent for only CGC genes
    :param df: df with raw pvalues
    :param path_weights: path to the weights
    :return:
    '''
    weight_dict = parse_optimized_weights(path_weights)
    weights = np.abs(np.array([weight_dict[m] for m in methods]))
    func = lambda x: trimmed_stouffer_w(x, weights) # lambda function to trimm p-values

    # Get stouffer weighted pvalues
    df['PVALUE_' + 'stouffer_w'] = df[['PVALUE_' + m for m in methods]].apply(func, axis=1)

    # Filter out nan p-values for fdr correction
    df = df[np.isfinite(df['PVALUE_' + 'stouffer_w'])].copy()
    df_nan = df[~np.isfinite(df['PVALUE_' + 'stouffer_w'])].copy()
    if df[np.isfinite(df['PVALUE_' + 'stouffer_w'])].shape[0] >0:
        df['QVALUE_' + 'stouffer_w'] = multipletests(df['PVALUE_' + 'stouffer_w'].values, method='fdr_bh')[1]
    else:
        df['QVALUE_' + 'stouffer_w'] = np.nan

    # Custom pvalue combination: Empirical Brown's Method -- including truncated method
    df = custom_combination(df, 'brown')

    # Custom pvalue combination: Fisher's Method -- including truncated method
    df = custom_combination(df, 'fisher')

    # Perform CGC correction
    cgc_set = load_cgc()
    df_cgc = df[df["SYMBOL"].isin(cgc_set)].copy()
    pvalues_cgc = df_cgc['PVALUE_' + 'stouffer_w'].values
    qvalues_cgc = multipletests(pvalues_cgc, method='fdr_bh')[1]
    df_cgc['QVALUE_CGC_' + 'stouffer_w'] = qvalues_cgc

    # Merge df with df_cgc dataframes
    df_final = pd.merge(left=df, right=df_cgc[["SYMBOL", "QVALUE_CGC_stouffer_w"]],
                        left_on="SYMBOL", right_on=["SYMBOL"], how="left")

    # Concat with non-corrected nan-containing dataframe
    df_final = pd.concat([df_final, df_nan],sort=True)
    return df_final


def partial_correction(df, fml_data):
    '''
    Perform a multipletest correction of the p-values discarding those nan p-values from oncodrivefml
    :param df:
    :param fml_data:
    :return:
    '''
    dh = pd.merge(left=df, right=fml_data[['SYMBOL', 'Q_VALUE',"SAMPLES","MUTS","MUTS_RECURRENCE"]], left_on=['SYMBOL'], right_on=['SYMBOL'], how="left")
    c = dh['Q_VALUE'].values
    mask = ~np.isnan(c)
    a = dh['PVALUE_' + 'stouffer_w'].values

    if len(a[mask]) == 0:
        print ("Warning: No data after filtering NaN OncodriveFML q-values")
        dh['QVALUE_' + 'stouffer_w'] = np.nan
        del dh['Q_VALUE']
        return dh
    c[mask] = multipletests(a[mask], method='fdr_bh')[1] # perform multiple test for those genes with at least two mutated samples
    dh['QVALUE_' + 'stouffer_w'] = c
    return dh


def include_excess(df,path_dndscv):
    '''

    :param df:
    :param path_dndscv:
    :return:
    '''
    dnds_data = pd.read_csv(path_dndscv, sep='\t', compression="gzip")
    columns = ['gene_name', 'wmis_cv',"wnon_cv","wspl_cv",'n_mis','n_non']
    if "wind_cv" in dnds_data.columns.values:
        columns = columns + ["wind_cv"]
    dh = pd.merge(left=df, right=dnds_data[columns], left_on=['SYMBOL'], right_on=['gene_name'], how="left")
    return dh


def filter_out_lowly_mutated(df,path_fml):
    '''
    Perform FDR correction of p-values for genes with at least two mutated samples.
    :param df: input dataframe with the p-values
    :param path_fml: path of FML, include information of mutated samples
    :return: input dataframe with a corrected p-value for samples with at least two mutated samples. nan otherwise.
    '''

    fml_data = pd.read_csv(path_fml, sep='\t', compression="gzip")
    dh = partial_correction(df, fml_data)
    return dh


def select_significant_bidders(row,QVALUE_threshold=0.1):
    '''
    Function to detect the significant bidders based on the
    :param row:
    :return:
    '''
    methods = []
    if len(str(row["All_Bidders"])) > 3: # Check whether there is a method (not nan)
        bidders = str(row["All_Bidders"]).split(",")
        for bidder in bidders:

            method_name = bidder.split("_")[0]
            name_key = "QVALUE_"+method_name
            if row[name_key] <= QVALUE_threshold:
                methods.append(method_name)
        return ",".join(methods)
    return ""


def add_significant_bidders(df):
    '''
    Add a column of the methods that have a significant bid
    :param df:
    :return: the input dataframe with a new column of significant bidders
    '''
    df["Significant_Bidders"] = df.apply(lambda row: select_significant_bidders(row),axis=1)
    return df


@click.command()
@click.option('--input_path', type=click.Path(exists=True), help="Path to input dataframe with SYMBOL and p-value for each method", required=True)
@click.option('--output_path', type=click.Path(), help="Path to output dataframe", required=True)
@click.option('--path_rankings', type=click.Path(), help="Path to dataframe produced by the voting system", required=True)
@click.option('--path_weights', type=click.Path(), help="Path to dataframe with weights", required=True)
@click.option('--path_fml', type=click.Path(), help="Path to OncodriveFML results folder", required=True)
@click.option('--path_dndscv', type=click.Path(), help="Path to dndsCV results folder", required=True)
@click.option('--brown', is_flag=True, help="append Empirical Brown combination")
@click.option('--fisher', is_flag=True, help="append Fisher combination")
def run_stouffer_script(input_path, output_path, path_rankings, path_weights, path_fml, path_dndscv, brown, fisher):

    # Read data
    df = pd.read_csv(input_path, sep='\t', compression="gzip")

    # Map with the ranking
    dg = retrieve_ranking(df, path_rankings)

    # Combine the pvalue
    dh = combine_pvals(dg, path_weights, METHODS)

    # Include the excess
    di = include_excess(dh, path_dndscv)

    # Add significant bidders
    dj = add_significant_bidders(di)

    # Remove q-values for samples with less than two mutated samples
    dk = filter_out_lowly_mutated(dj, path_fml)

    # Display the table with sorted columns
    pqvals = reduce(lambda x, y: x + y, map(lambda x: ['PVALUE_' + x, 'QVALUE_' + x], METHODS))
    column_order = ["SYMBOL"] + pqvals + ["PVALUE_stouffer_w", "QVALUE_stouffer_w", "QVALUE_CGC_stouffer_w"]
    if brown:
        column_order += list(map(lambda x: x + '_brown', ['PVALUE', 'QVALUE', 'PVALUE_trunc', 'QVALUE_trunc']))
    if fisher:
        column_order += list(map(lambda x: x + '_fisher', ['PVALUE', 'QVALUE', 'PVALUE_trunc', 'QVALUE_trunc']))
    column_order += ["All_Bidders", "Significant_Bidders", "Median_Ranking", "RANKING", "Total_Bidders"] + \
                    ["wmis_cv", "wnon_cv", "wspl_cv", "SAMPLES", "MUTS", "MUTS_RECURRENCE", "n_mis", "n_non"]

    if "wind_cv" in dg.columns.values:
        column_order.append("wind_cv")

    if "RANKING_BORDA" in dg.columns.values:
        column_order.append("RANKING_BORDA")

    dk[column_order].to_csv(output_path, sep='\t', index=False, compression="gzip")


if __name__ == '__main__':
    run_stouffer_script()
