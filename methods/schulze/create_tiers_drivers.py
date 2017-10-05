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

@click.command()
@click.option('--input_dir',type=click.Path(exists=True),help="Directory of the dir to be parsed",required=True)
@click.option('--output_file', type=click.Path(),required=True)
@click.option('--threshold', help="Directory of the output reports",required=False,default=0.01,type=float)
@click.option('--name_method', help="Name of the method in the output dataframe",default="SCHULZE_THRESHOLD_STOUFFER_WEIGHTED")
@click.option('--column_filter', help="Column to be used by the filtering",default="QVALUE_schulze_weighted")
@click.option('--pattern', help="Pattern look in the input_dir",default="_combined.tsv")


def run_create_tiers(input_dir,output_file,threshold,name_method,column_filter,pattern):
    l_cancers = []


    for filename in glob.iglob(input_dir+"/*"+pattern):
        print ("Working on:", filename)
        m = re.search('([A-Z]+)'+pattern, filename)
        if m:
            cancer_type = m.group(1)
            print (cancer_type)
            df = pd.read_csv(filename,sep="\t")
            df.sort_values(column_filter,inplace=True)
            df_f = df[~np.isnan(df[column_filter])&(df[column_filter]<0.5)].copy()
            ranking_limit = df_f.sort_values("RANKING",ascending=False).head(1)["RANKING"].values[0]
            if (ranking_limit):
                dfq = df_f[df_f["RANKING"]<ranking_limit].copy()
                dfq = dfq[np.isfinite(dfq[column_filter])].copy()
                d_class_3tiers = classify_genes_tiers(dfq,column_filter=column_filter,threshold=threshold)
                dfq["TIER"] = dfq.apply(lambda row: d_class_3tiers[row["SYMBOL"]],axis=1)
                #dfq["TIER"] = dfq.apply(lambda row: "3" if (row["TIER"]==2 and not("," in row["All_Bidders"])) else row["TIER"],axis=1)
                dfq["METHOD_NAME"] = name_method
                l_cancers.append(dfq[["SYMBOL","Cancer_Type","METHOD_NAME","TIER","All_Bidders",column_filter,"RANKING"]])
    df_tiers = pd.concat(l_cancers)
    df_tiers.to_csv(output_file,sep="\t",index=False)
    return df_tiers




if __name__ == '__main__':


    run_create_tiers()



