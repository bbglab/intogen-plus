import os

import click as click
import pandas as pd
import numpy as np
from scipy.stats import combine_pvalues, uniform
from statsmodels.stats.multitest import multipletests


#DEFAULT_METHODS = ['oncodrivefml', 'oncodriveclust', 'oncodriveomega', 'hotmapssignature', 'mutsigcv']
DEFAULT_METHODS = ["oncodriveclustl", "dndscv","oncodrivefml", "hotmapssignature","edriver","cbase"]
column_keys = {
        "hotmapssignature":     ["GENE",        "q-value",      "Min p-value"],
        "oncodrivefml":         ["SYMBOL",      "Q_VALUE",      "P_VALUE"],
        "dndscv":               ["gene_name",   "qallsubs_cv",  "pallsubs_cv"],
        "edriver":              ["SYMBOL",      "QVALUE",       "PVALUE"],
        "cbase":                ["gene",        "q_phi_pos",    "p_phi_pos"],
        "oncodriveclustl":      ["SYM",         "E_QVAL",       "E_PVAL"]
        # "oncodriveclust":       ["SYMBOL",      "QVALUE",       "PVALUE"],
        # "mutsigcv": ["gene","q", "p"],
        # "oncodriveomega": ["SYMBOL","q_value", "p_value"],
    }
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

    nan_mask = np.isnan(pvals) # used to exclude nan p-values from combination
    anti_one_mask = pvals[pvals < 1.] # used to exclude p-values 1.0 from combination
    mask = ~nan_mask & anti_one_mask
    return pvals[mask], mask


def truncate(pvals, threshold=1e-16):

    mask = (pvals < threshold)
    pvals[mask] = threshold
    return pvals


def trimmed_stouffer_w(pvals, weights):
    """
    conducts stouffer_w where pvals and weights are clean from nans
    """

    reduced_pvals, mask = trim_nans(pvals)
    reduced_weights = weights[mask]
    return stouffer_w(truncate(reduced_pvals), weights=reduced_weights)


def load_cgc():
    '''
    Loads the CGC set and returns a set of CGC genes
    :return: set of cgc gene names
    '''
    df_cgc = pd.read_csv(os.path.join(os.environ['SCHULZE_DATA'], "CGC_set.tsv"), sep="\t")
    return set(df_cgc["Gene Symbol"].values)


def combine_pvals(df, path_weights):
    '''

    :param df: df with raw pvalues
    :param path_weights: path to the weights
    :return:
    '''

    weight_dict = parse_optimized_weights(path_weights)
    weights = np.abs(np.array([weight_dict[m] for m in DEFAULT_METHODS]))

    func = lambda x: trimmed_stouffer_w(x, weights)
    # Get the stouffer pvalues
    df['PVALUE_' + 'stouffer_w'] = df[['PVALUE_' + m for m in DEFAULT_METHODS]].apply(func, axis=1)
    # Filter out nan p-values for fdr correction
    df_non_nan = df[np.isfinite(df['PVALUE_' + 'stouffer_w'])].copy()
    df_non_nan['QVALUE_' + 'stouffer_w'] = multipletests(df_non_nan['PVALUE_' + 'stouffer_w'].values, method='fdr_bh')[1]
    # Perform cgc correction
    cgc_set = load_cgc()
    df_cgc = df_non_nan[df_non_nan["SYMBOL"].isin(cgc_set)].copy()
    pvalues_cgc = df_cgc['PVALUE_' + 'stouffer_w'].values
    qvalues_cgc = multipletests(pvalues_cgc, method='fdr_bh')[1]
    df_cgc['QVALUE_CGC_' + 'stouffer_w'] = qvalues_cgc
    # Merge with the non_nan dataframe
    df_final_non_nan = pd.merge(left=df_non_nan,right=df_cgc[["SYMBOL","QVALUE_CGC_stouffer_w"]],left_on="SYMBOL",right_on=["SYMBOL"],how="left")
    # Concat the non_nan dataframe with the non-corrected nan-containing dataframe
    df_final_nan = df[~np.isfinite(df['PVALUE_' + 'stouffer_w'])].copy()
    df_final = pd.concat([df_final_non_nan,df_final_nan])
    return df_final


def partial_correction(df, fml_data):
    '''

    :param df:
    :param fml_data:
    :return:
    '''
    dh = pd.merge(left=df, right=fml_data[['SYMBOL', 'Q_VALUE',"SAMPLES","MUTS","MUTS_RECURRENCE"]], left_on=['SYMBOL'], right_on=['SYMBOL'], how="left")
    c = dh['Q_VALUE'].values
    mask = ~np.isnan(c)
    a = dh['PVALUE_' + 'stouffer_w'].values

    if len(a[mask]) == 0:
        raise RuntimeError("No data after filtering NaN OncodriveFML q-values")

    c[mask] = multipletests(a[mask], method='fdr_bh')[1]
    dh['QVALUE_' + 'stouffer_w'] = c
    del dh['Q_VALUE']
    return dh

def include_excess(df,path_dndscv):
    '''

    :param df:
    :param path_dndscv:
    :return:
    '''
    dnds_data = pd.read_csv(path_dndscv, sep='\t', compression="gzip")
    columns = ['gene_name', 'wmis_cv',"wnon_cv","wspl_cv"]
    if "wind_cv" in dnds_data.columns.values:
        columns = columns + ["wind_cv"]
    dh = pd.merge(left=df, right=dnds_data[columns], left_on=['SYMBOL'], right_on=['gene_name'], how="left")
    return dh

def combine_from_tumor(df, path_to_output, path_fml):
    '''

    :param df:
    :param path_to_output:
    :param path_fml:
    :return:
    '''

    fml_data = pd.read_csv(path_fml, sep='\t', compression="gzip")
    dh = partial_correction(df, fml_data)
    column_order = ["SYMBOL","PVALUE_cbase","PVALUE_dndscv","QVALUE_dndscv","PVALUE_edriver","QVALUE_edriver","PVALUE_hotmapssignature","QVALUE_hotmapssignature","PVALUE_oncodriveclustl","QVALUE_oncodriveclustl","PVALUE_oncodrivefml","QVALUE_oncodrivefml","PVALUE_stouffer_w","QVALUE_stouffer_w","QVALUE_CGC_stouffer_w","All_Bidders","Significant_Bidders","Median_Ranking","RANKING","Total_Bidders","wmis_cv","wnon_cv","wspl_cv","SAMPLES","MUTS","MUTS_RECURRENCE"]
    if "wind_cv" in dh.columns.values:
        column_order.append("wind_cv")
    dh[column_order].to_csv(path_to_output, sep='\t', index=False, compression="gzip")


def select_significant_bidders(row,QVALUE_threshold=0.1):
    '''
    Function to detect the significant bidders based on the
    :param row:
    :return:
    '''
    methods = []
    if len(str(row["All_Bidders"]))>3: # Check whether there is a method (not nan)
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

def run_stouffer_script(input_path, output_path, path_rankings, path_weights, path_fml,path_dndscv):

    df = pd.read_csv(input_path, sep='\t', compression="gzip")
    dg = retrieve_ranking(df, path_rankings)
    dh = combine_pvals(dg, path_weights)
    di = include_excess(dh, path_dndscv)
    dg = add_significant_bidders(di)
    combine_from_tumor(dg, output_path, path_fml)


if __name__ == '__main__':
    run_stouffer_script()

