'''
Prepare vetting files to be used for the calling of drivers script
'''
import click
import os
from collections import defaultdict
import glob
import re
import pandas as pd
import numpy as np


artifacts=['HMCN1', 'TTN', 'OBSCN', 'GPR98', "RYR2","RYR3"]

def count_percentage_signature(grp):
    signatures_hypermutators = ["Signature.9", "Signature.10"]
    total = len(grp)
    d = defaultdict(int)
    d_norm = defaultdict(float)
    for sig in grp:
        d[sig] += 1
    for s in d.keys():
        d_norm[s] = (d[s] / total)
    return d_norm[signatures_hypermutators[0]], d_norm[signatures_hypermutators[1]]


def assign_type_mut(value):
    # C[C>T]G SNP; T[C>-]A INDEL; T[C>CA]A INDEL
    muts = value[2:-2]
    muts = muts.replace("-", "")
    nts = muts.split(">")
    if len(nts[0]) != len(nts[1]):
        return "INDEL"
    else:
        return "SNP"


def get_percentage_palindrome(grp):
    n_pals = 0
    n_neutral = 0
    for e in grp:
        if e == "yes":
            n_pals += 1
        else:
            n_neutral += 1
    return n_pals / (n_pals + n_neutral)

def get_warning_muts_sample(df,n_muts_gene):
    '''

    :param df:
    :return:
    '''
    muts_sample_gene = df.groupby(["SAMPLE", "GENE"], as_index=False).agg({"POSITION": "count"})
    genes_warning_samples = muts_sample_gene[muts_sample_gene["POSITION"] >= n_muts_gene].groupby("GENE",as_index=False).agg({"SAMPLE": "count"})
    genes_warning_samples.rename(columns={"SAMPLE": "Samples_3muts"}, inplace=True)
    return genes_warning_samples


def analysis_signatures_gene(file_sigs,df,cohort):
    '''

    :param file_sigs:
    :return:
    '''
    if os.path.exists(file_sigs):
        df_sigs = pd.read_csv(file_sigs, sep="\t")
        df_sigs.fillna(0.0, inplace=True)
        # get the signature with the highest contribution for a particular sample
        df_sigs["signature_max"] = df_sigs.apply(
            lambda row: row.index[list(row.values).index(np.nanmax(row.values[2:]))], axis=1)
        df_combined = pd.merge(df, df_sigs[["Mutation_type", "signature_max", "Sample"]],
                               left_on=["SAMPLE", "MUTATION_TYPE"], right_on=["Sample", "Mutation_type"],
                               how="left")
        df_combined["signature_max"].fillna("", inplace=True)
        df_combined.drop(columns=["Mutation_type", "Sample"], inplace=True)
    else:
        df_combined = df.copy()
        df_combined["signature_max"] = ""
    df_combined["TYPE_MUT"] = df_combined.apply(lambda row: assign_type_mut(row["MUTATION_TYPE"]), axis=1)
    # Assign % of signatures per gene

    df_genes = df_combined[df_combined["TYPE_MUT"] == "SNP"].groupby("GENE", as_index=False).agg(
        {"signature_max": count_percentage_signature, "POSITION": "count"}).sort_values("POSITION",
                                                                                        ascending=False)
    df_genes["COHORT"] = cohort
    df_genes["Signature9"] = df_genes.apply(lambda row: row["signature_max"][0], axis=1)
    df_genes["Signature10"] = df_genes.apply(lambda row: row["signature_max"][1], axis=1)
    df_genes.drop(columns=["POSITION", "signature_max"], inplace=True)
    return df_genes,df_combined


def get_ratio_indels(df_combined):
    '''

    :param df_:
    :return:
    '''
    df_numbers = df_combined.groupby(["GENE", "TYPE_MUT"], as_index=False).agg({"POSITION": "count"})
    df_agg = df_numbers.pivot_table(columns=["TYPE_MUT"], values=["POSITION"], index=["GENE"],
                                    fill_value=0.0).reset_index()
    df_agg.columns = df_agg.columns.droplevel()
    if df_numbers[df_numbers["TYPE_MUT"] == "INDEL"].shape[0] == 0:
        df_agg["INDEL"] = 0.0
        df_agg.columns = ["GENE", "SNP", "INDEL"]
    else:

        df_agg.columns = ["GENE", "INDEL", "SNP"]

    df_agg["INDEL/SNP"] = df_agg.apply(lambda row: row["INDEL"] / (row["SNP"] + row["INDEL"]), axis=1)

    return df_agg

def check_expression(row,d_expression_tcga):

    if row["cancer_type"] in d_expression_tcga:

        return (row["GENE"] in d_expression_tcga[row["cancer_type"]])
    else:
        return (row["GENE"] in d_expression_tcga["PANCANCER"])


def filter_by_expression(df, expresison_file_pcawg, expresison_file_tcga):
    '''
    Filter dataframe df by expression of genes
    :param df:
    :return:
    '''
    # Load expression from TCGA
    df_expression_tcga = pd.read_csv(expresison_file_tcga, sep="\t", names=["Cancer_Type", "GENES"])
    d_expression_tcga = {}
    for index, row in df_expression_tcga.drop_duplicates().iterrows():
        d_expression_tcga[row["Cancer_Type"]] = row["GENES"].split(",")

    df["Warning_Expression"] = df.apply(lambda row: check_expression(row,d_expression_tcga),axis=1)
    return df

def filter_by_polymorphism(df,file_exact):
    '''
    Filter by population polymorphism
    :param df:
    :return:
    '''
    df_exac = pd.read_csv(file_exact,sep="\t")
    df_exac = df_exac[df_exac["canonical"]]
    df_exac_filtered = df_exac[["gene", "oe_syn", "oe_syn", "oe_mis", "bp"]].drop_duplicates()
    df_final_total = pd.merge(df_exac_filtered, df, left_on=["gene"], right_on=["GENE"],how="right")
    df_final_total.drop(columns=["gene"], inplace=True)
    df_final_total["Warning_Germline"] = df_final_total.apply(lambda row: row["oe_syn"] > 1.5 or row["oe_mis"] > 1.5 or row["oe_lof"] > 1.5, axis=1)
    return df_final_total

def filter_by_gene_lenght(df,zscore_threshold=5):
    '''

    :param df:
    :param zscore_threshold:
    :return:
    '''
    mean = np.nanmean(df["bp"].values)
    std = np.nanstd(df["bp"].values)
    df["Z-score length"] = df.apply(lambda row: (row["bp"] - mean) / std, axis=1)
    df["Warning_size"] = df.apply(lambda row: row["Z-score length"] > zscore_threshold, axis=1)
    return df

def filter_by_olfactory_receptors(df, of_file):
    '''
    Check whether is an olfactory receptor
    :param df:
    :param of_file:
    :return:
    '''
    df_blacklisted = pd.read_csv(of_file, sep="\t")
    orfs = list(df_blacklisted["Symbol"].unique())
    df["OR_Warning"] = df.apply(lambda row: True if row["GENE"] in orfs else False,axis=1)
    return df

def main(paths,pattern,n_muts_gene,info_cohorts,output,expression_file_tcga,exact_germline_mutations,of_file):

    l = []
    for input_path in paths:
        for filein in glob.glob(input_path + "*.in.gz"):
            df = pd.read_csv(filein, sep="\t")
            m = re.search(pattern, filein)
            cohort = m.group(1)
            if "PANEL" in cohort: # This is not intended for panels is just to vett large projects
                continue
            # 1. Samples with more than 2 mutations in a gene are likely an artifact
            genes_warning_samples = get_warning_muts_sample(df,n_muts_gene)
            # 2. Analysis of signatures
            file_sigs = os.path.join(input_path, cohort + ".out.gz", "mutation_sign_prob.tsv")
            df_genes,df_combined=analysis_signatures_gene(file_sigs,df,cohort)
            # 3. Ratio Ratio indels SNP
            df_agg = get_ratio_indels(df_combined)
            df_agg["COHORT"] = cohort
            # 4. Combine all information
            df_info_total = pd.merge(df_agg, df_genes, how="outer")
            df_info_total = pd.merge(df_info_total, genes_warning_samples, how="left")
            df_info_total.fillna(0.0, inplace=True)
            l.append(df_info_total.copy())



    # Combine all data into a dataframe
    df_final_signatures = pd.concat(l,sort=True)
    # Get info cohorts
    df_info_cohort = pd.read_csv(info_cohorts, sep="\t")
    # Add information of cancer type
    df_final_signatures_info = pd.merge(df_final_signatures, df_info_cohort)
    # Filter by expression
    df_final_signatures_info=filter_by_expression(df_final_signatures_info,expression_file_tcga)
    # Filter by Polymorphism
    df_final_signatures_info=filter_by_polymorphism(df_final_signatures_info,exact_germline_mutations)
    # Add filter by Gene Lenght
    df_final_signatures_info = filter_by_gene_lenght(df_final_signatures_info)
    # Add filter by olfactory receptors
    df_final_signatures_info = filter_by_olfactory_receptors(df_final_signatures_info,of_file)
    # Add filter by known artifacts
    df_final_signatures_info["Warning_Artifact"] = df_final_signatures_info.apply(lambda row: row["GENE"] in artifacts,axis=1)
    # Save it
    df_final_signatures_info.to_csv(output,sep="\t",index=False,compression="gzip")



@click.command()
@click.option('-i', '--intogen', 'intogen',type=click.Path(),  help='Path to intogen data',required=True)
@click.option('-w', '--hartwig', 'hartwig', type=click.Path(),  help='Path to hartwig data', required=True)
@click.option('-s', '--stjude', 'stjude', type=click.Path(),  help="Path to stjude data", required=True)
@click.option('-c', '--info_cohorts', 'info_cohorts', type=click.Path(),  help="Path to info of the cohorts", required=True)
@click.option('-o', '--output', 'output', type=click.Path(),  help="Path to the output. Formatted in .tsv.gz", required=True)
@click.option('-t', '--expression_file_tcga', 'exp_tcga', type=click.Path(),  help="Path to non expressed genes form tcga", required=True)
@click.option('-e', '--exact_germline_mutations', 'exact', type=click.Path(),  help="Path to exact germline variants", required=True)
@click.option('-r', '--olfatory_receptors', 'olfatory_receptors', type=click.Path(),  help="Path to the olfactory receptors", required=True)
def cmdline(intogen, hartwig, stjude,info_cohorts,output,exp_tcga,exact, olfatory_receptors):
    paths = [intogen,hartwig,stjude]
    # configuration
    pattern = re.compile("\/([A-Z0-9_\-\,]+)\.in\.gz")
    n_muts_gene = 3
    main(paths,pattern,n_muts_gene,info_cohorts,output,exp_tcga,exact,olfatory_receptors)


if __name__ == "__main__":
    cmdline()






