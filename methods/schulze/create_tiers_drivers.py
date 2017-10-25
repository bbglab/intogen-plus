import csv
import gzip

import pandas as pd
import glob
import re
import numpy as np
import click

def classify_genes_tiers(df,threshold=0.01,column_ranking="RANKING",column_filter = "QVALUE_stouffer_weighted"):
    df.sort_values("RANKING",ascending=False,inplace=True)
    found = False
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
            d_class[row["SYMBOL"]] = 2
        elif row["RANKING"] < threshold_rejected and row[column_filter] > threshold:
            d_class[row["SYMBOL"]] = 4
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
    l = set(df[df[column]<threshold]["SYMBOL"].values)
    return l
def rescue_genes(row,list_genes_recovered):
    if row["TIER"] <4:
        return row["TIER"]
    if row["SYMBOL"] in list_genes_recovered:
        return 3
    return 4



@click.command()
@click.option('--input',type=click.Path(exists=True),help="File to be parsed",required=True)
@click.option('--output_file', type=click.Path(),required=True)
@click.option('--threshold', help="Directory of the output reports",required=False,default=0.01,type=float)
@click.option('--name_method', help="Name of the method in the output dataframe",default="SCHULZE_THRESHOLD_STOUFFER_WEIGHTED")
@click.option('--column_filter', help="Column to be used by the filtering",default="QVALUE_schulze_weighted")
@click.option('--column_filter_cgc', help="Column to be used by the rescue of CGC",default="QVALUE_CGC_stouffer_w")


def run_create_tiers(input, output_file, threshold, name_method, column_filter,column_filter_cgc):

    df = pd.read_csv(input, sep="\t", compression="gzip")
    df.sort_values(column_filter,inplace=True)
    df_f = df[~np.isnan(df[column_filter])&(df[column_filter]<0.5)].copy()

    ranking_limit = df_f.sort_values("RANKING",ascending=False).head(1)["RANKING"].values[0] if len(df_f) > 1 else None

    headers = ["SYMBOL", "METHOD_NAME","TIER","All_Bidders", column_filter, "RANKING"]
    if ranking_limit:

        # Compute tiers
        dfq = df_f[df_f["RANKING"] < ranking_limit].copy()
        dfq = dfq[np.isfinite(dfq[column_filter])].copy()
        d_class_3tiers = classify_genes_tiers(dfq,column_filter=column_filter,threshold=threshold)
        dfq["TIER"] = dfq.apply(lambda row: d_class_3tiers[row["SYMBOL"]],axis=1)
        dfq["METHOD_NAME"] = name_method
        rescued_genes = get_recovered_genes(dfq,column_filter_cgc,threshold)
        dfq["TIER"] = dfq.apply(lambda row: rescue_genes(row,rescued_genes),axis=1)
        df_tiers = dfq[headers]
        df_tiers.to_csv(output_file, sep="\t", index=False, compression="gzip")
    else:

        # No results
        with gzip.open(output_file, 'wt') as fd:
            csv.writer(fd, delimiter='\t').writerow(headers)


if __name__ == '__main__':

    run_create_tiers()



