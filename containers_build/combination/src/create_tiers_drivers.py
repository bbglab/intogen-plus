import csv
import gzip
import math

import click
import numpy as np
import pandas as pd


def classify_genes_tiers(df, threshold=0.01, column_filter="QVALUE_stouffer_w"):
    '''
    Split the input dataframe into several classes according to their ranking and their p-value
    :param df: input dataframe to be classified
    :param threshold: threshold for the tier 1
    :param column_filter: name of the column to be filtered
    :return: a dictionary whose keys are gene symbols and values are the tiers
    '''

    df.sort_values("RANKING",ascending=False,inplace=True)

    threshold_rejected = df.shape[0]
    threshold_accepted = 1
    d_class = {}

    for index, row in df.iterrows():
        if row[column_filter] <= threshold:
            threshold_rejected = row["RANKING"] + 1
            break
    df.sort_values("RANKING",ascending=True,inplace=True)
    for index, row in df.iterrows():
        if row[column_filter] > threshold:
            threshold_accepted= row["RANKING"] - 1
            break

    for index, row in df.iterrows():
        if row["RANKING"] <= threshold_accepted:
            d_class[row["SYMBOL"]] = 1
        elif row["RANKING"] >= threshold_rejected:
            d_class[row["SYMBOL"]] = 4
        elif row["RANKING"] < threshold_rejected and row[column_filter] <= threshold:
            d_class[row["SYMBOL"]] = 3
        else:
            d_class[row["SYMBOL"]] = 4
    return d_class


def get_recovered_genes(df,column,threshold):
    '''
    Method to rescue CGC genes after FDR correction limited to CGC genes.
    :param df: the input dataframe
    :param column: colum of the QVALUE of FDR-CGC correction
    :param threshold: the limit threshold
    :return: the symbol ID of the genes below the threshold
    '''
    list_genes_recovered = set(df[df[column]<threshold]["SYMBOL"].values)
    return list_genes_recovered


def rescue_genes(row,list_genes_recovered):
    '''

    :param row:
    :param list_genes_recovered:
    :return:
    '''
    if row["TIER"] == 1:
        return row["TIER"]
    if row["SYMBOL"] in list_genes_recovered:
        return 2
    return row["TIER"]




def set_role(data, distance_threshold=0.1):
    """Set the role according to the DNDS output"""
    if data['wmis_cv'] < 1 and data['wnon_cv'] < 1:  # threshold
        return "ambiguous"
    # Check wmis
    wmis = data['wmis_cv']
    if wmis >= 1 and data["n_mis"] == 0:
        wmis = 1

    # Check wnon
    wnon = data['wnon_cv']
    if wnon >= 1 and data["n_non"] == 0:
        wnon = 1
    # Those cases with w_non and w_mis <=1 are not informative
    if wnon <= 1 and wmis <= 1:
        return "ambiguous"

    distance = (wmis - wnon) / math.sqrt(2)
    if distance_threshold is not None and abs(distance) < distance_threshold:
        return "ambiguous"
    else:
        if distance > 0:
            return 'Act'
        elif distance < 0:
            return 'LoF'
        else:
            return "ambiguous"

@click.command()
@click.option('--input',type=click.Path(exists=True),help="File to be parsed",required=True)
@click.option('--output_file', type=click.Path(),required=True)
@click.option('--threshold', help="Threshold to be used to filter all genes",required=False,default=0.01,type=float)
@click.option('--threshold_cgc', help="Threshold to be used to filter CGC genes",required=False,default=0.25,type=float)
@click.option('--column_filter', help="Column to be used by the filtering",default="QVALUE_stouffer_w")
@click.option('--column_filter_cgc', help="Column to be used by the rescue of Cancer Gene Census genes",default="QVALUE_CGC_stouffer_w")
def run_create_tiers(input, output_file, threshold, threshold_cgc, column_filter,column_filter_cgc):

    df = pd.read_csv(input, sep="\t", compression="gzip")
    df.sort_values(column_filter,inplace=True)
    df_f = df[~np.isnan(df[column_filter])&(df[column_filter]<0.5)].copy() # Select only a portion of likely candidates, make the ranking faster
    ranking_limit = df_f.sort_values("RANKING",ascending=False).head(1)["RANKING"].values[0] if len(df_f) > 1 else None
    headers = ["SYMBOL", "TIER","All_Bidders","Significant_Bidders", column_filter, "RANKING","MUTS","SAMPLES",'wmis_cv',"wnon_cv","wspl_cv","n_mis","n_non"] # ,

    if ranking_limit:
        # Compute tiers
        dfq = df_f[df_f["RANKING"] < ranking_limit].copy() # select those positions before the limit
        dfq = dfq[np.isfinite(dfq[column_filter])].copy() # with finite q-value
        d_class_3tiers = classify_genes_tiers(dfq,column_filter=column_filter,threshold=threshold) # perform the classification
        dfq["TIER"] = dfq.apply(lambda row: d_class_3tiers[row["SYMBOL"]],axis=1)
        rescued_genes = get_recovered_genes(dfq,column_filter_cgc,threshold_cgc) # perform the rescue of cgc genes
        dfq["TIER"] = dfq.apply(lambda row: rescue_genes(row,rescued_genes),axis=1)
        df_tiers = dfq[headers]
        df_tiers['ROLE'] = df_tiers.apply(lambda row: set_role(row), axis=1)
        df_tiers.to_csv(output_file, sep="\t", index=False, compression="gzip")
    else:
        # No results
        with gzip.open(output_file, 'wt') as fd:
            csv.writer(fd, delimiter='\t').writerow(headers)


if __name__ == '__main__':
    run_create_tiers()



