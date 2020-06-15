import os

import pandas as pd


def filter_samples_by_nmuts(df, n_muts_gene):
    muts_sample_gene = df.groupby(["SAMPLE", "GENE"], as_index=False).agg({"POSITION": "count"})
    genes_warning_samples = muts_sample_gene[muts_sample_gene["POSITION"] >= n_muts_gene].groupby("GENE", as_index=False).agg({"SAMPLE": "count"})
    genes_warning_samples.rename(columns={"SAMPLE": f"Samples_{n_muts_gene}muts"}, inplace=True)
    return genes_warning_samples


# TODO ensure our ctypes match TCGA acronyms
def check_expression(row, ctype, d_expression_tcga):
    if ctype in d_expression_tcga:
        return row["GENE"] in d_expression_tcga[ctype]
    else:
        return row["GENE"] in d_expression_tcga["PANCANCER"]


def filter_by_expression(df, ctype):
    """Filter dataframe df by expression of genes"""
    # Load expression from TCGA
    expresison_file_tcga = os.path.join(os.environ['INTOGEN_DATASETS'], 'others', 'non_expressed_genes_tcga.tsv')
    df_expression_tcga = pd.read_csv(expresison_file_tcga, sep="\t", names=["Cancer_Type", "GENES"])
    d_expression_tcga = {}
    for index, row in df_expression_tcga.drop_duplicates().iterrows():
        d_expression_tcga[row["Cancer_Type"]] = row["GENES"].split(",")
    df["Warning_Expression"] = df.apply(lambda row: check_expression(row, ctype, d_expression_tcga), axis=1)
    return df


def filter_by_polymorphism(df):
    """Filter by population polymorphism"""
    file_exact = os.path.join(os.environ['INTOGEN_DATASETS'], 'postprocess', 'constraint.txt.gz')
    df_exac = pd.read_csv(file_exact,sep="\t")
    df_exac = df_exac[df_exac["canonical"]]
    df_exac_filtered = df_exac[["gene", "oe_syn", "oe_lof", "oe_mis"]].drop_duplicates()
    df_final_total = pd.merge(df_exac_filtered, df, left_on=["gene"], right_on=["GENE"], how="right")
    df_final_total.drop(columns=["gene"], inplace=True)
    df_final_total[["oe_syn", "oe_mis", "oe_lof"]].fillna(0.0,inplace=True)
    df_final_total["Warning_Germline"] = df_final_total.apply(lambda row: row["oe_syn"] > 1.5 or row["oe_mis"] > 1.5 or row["oe_lof"] > 1.5, axis=1)
    return df_final_total


def filter_by_olfactory_receptors(df):
    """Check whether is an olfactory receptor"""
    of_file = os.path.join(os.environ['INTOGEN_DATASETS'], 'others', 'olfactory_receptors.tsv')
    df_blacklisted = pd.read_csv(of_file, sep="\t")
    orfs = list(df_blacklisted["Symbol"].unique())
    df["OR_Warning"] = df.apply(lambda row: True if row["GENE"] in orfs else False,axis=1)
    return df
