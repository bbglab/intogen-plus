
import glob
import json
import logging
import re
import sys
import tempfile
from collections import defaultdict
from os import path

import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None

FOLDER = path.dirname(path.abspath(__file__))


def count_percentage_signature(grp):
    signatures_hypermutators = ["Signature.9", "Signature.10"]
    total = len(grp)
    d = defaultdict(int)
    for sig in grp:
        d[sig] += 1
    d_norm = defaultdict(float)
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


def filter_samples_by_nmuts(df, n_muts_gene):
    muts_sample_gene = df.groupby(["SAMPLE", "GENE"], as_index=False).agg({"POSITION": "count"})
    genes_warning_samples = muts_sample_gene[muts_sample_gene["POSITION"] >= n_muts_gene].groupby("GENE", as_index=False).agg({"SAMPLE": "count"})
    genes_warning_samples.rename(columns={"SAMPLE": "Samples_3muts"}, inplace=True)
    return genes_warning_samples


def analysis_signatures_gene(file_sigs, df):
    if path.exists(file_sigs):
        df_sigs = pd.read_csv(file_sigs, sep="\t")
        df_sigs.fillna(0.0, inplace=True)
        if 'Unnamed: 0' in df_sigs.columns.values:
            df_sigs.drop(["Unnamed: 0"], axis=1, inplace=True)
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
    df_genes["Signature9"] = df_genes.apply(lambda row: row["signature_max"][0], axis=1)
    df_genes["Signature10"] = df_genes.apply(lambda row: row["signature_max"][1], axis=1)
    df_genes.drop(columns=["POSITION", "signature_max"], inplace=True)
    return df_genes, df_combined


def get_ratio_indels(df_combined):
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
    if row["CANCER_TYPE"] in d_expression_tcga:
        return row["GENE"] in d_expression_tcga[row["CANCER_TYPE"]]
    else:
        return row["GENE"] in d_expression_tcga["PANCANCER"]


def filter_by_expression(df, expresison_file_tcga):
    """Filter dataframe df by expression of genes"""
    # Load expression from TCGA
    df_expression_tcga = pd.read_csv(expresison_file_tcga, sep="\t", names=["Cancer_Type", "GENES"])
    d_expression_tcga = {}
    for index, row in df_expression_tcga.drop_duplicates().iterrows():
        d_expression_tcga[row["Cancer_Type"]] = row["GENES"].split(",")
    df["Warning_Expression"] = df.apply(lambda row: check_expression(row, d_expression_tcga), axis=1)
    return df


def filter_by_polymorphism(df, file_exact):
    """Filter by population polymorphism"""
    df_exac = pd.read_csv(file_exact,sep="\t")
    df_exac = df_exac[df_exac["canonical"]]
    df_exac_filtered = df_exac[["gene", "oe_syn", "oe_lof", "oe_mis"]].drop_duplicates()
    df_final_total = pd.merge(df_exac_filtered, df, left_on=["gene"], right_on=["GENE"],how="right")
    df_final_total.drop(columns=["gene"], inplace=True)
    df_final_total[["oe_syn","oe_mis","oe_lof"]].fillna(0.0,inplace=True)
    df_final_total["Warning_Germline"] = df_final_total.apply(lambda row: row["oe_syn"] > 1.5 or row["oe_mis"] > 1.5 or row["oe_lof"] > 1.5, axis=1)
    return df_final_total


def filter_by_olfactory_receptors(df, of_file):
    """Check whether is an olfactory receptor"""
    df_blacklisted = pd.read_csv(of_file, sep="\t")
    orfs = list(df_blacklisted["Symbol"].unique())
    df["OR_Warning"] = df.apply(lambda row: True if row["GENE"] in orfs else False,axis=1)
    return df


def include_literature(df, cancermine_file):
    # read
    cancermine = pd.read_csv(cancermine_file, sep="\t")
    cancermine_g = cancermine.groupby("gene_normalized", as_index=False).agg({"pmid": "count"})
    cancermine_g.rename(columns={"pmid": "n_papers", "gene_normalized": "GENE"}, inplace=True)
    # Match with genes
    df = df.merge(cancermine_g, how="left")
    df["n_papers"].fillna(0, inplace=True)
    return df


PATTERN = re.compile("\/([A-Z0-9_\-\,]+)\.in\.gz")


def prepare(paths, exp_tcga, olfatory_receptors, cancermine, germline, info_cohorts, muts=3):
    """Prepare vetting files to be used for the calling of drivers script"""

    l = []
    for path_ in paths:
        for filein in glob.glob(path.join(path_, "*.in.gz")):
            df = pd.read_csv(filein, sep="\t")
            m = PATTERN.search(filein)
            cohort = m.group(1)
            # 1. Samples with more than 2 mutations in a gene are likely an artifact
            genes_warning_samples = filter_samples_by_nmuts(df, muts)
            # 2. Analysis of signatures
            file_sigs = path.join(path_, cohort + ".signature_likelihood")
            df_genes, df_combined = analysis_signatures_gene(file_sigs, df)
            df_genes["COHORT"] = cohort
            # 3. Ratio Ratio indels SNP
            df_agg = get_ratio_indels(df_combined)
            df_agg["COHORT"] = cohort
            # 4. Combine all information
            df_info_total = pd.merge(df_agg, df_genes, how="outer")
            df_info_total = pd.merge(df_info_total, genes_warning_samples, how="left")
            df_info_total.fillna(0.0, inplace=True)
            l.append(df_info_total)

    # Combine all data into a dataframe
    df_final_signatures = pd.concat(l, sort=True)
    # Get info cohorts
    df_info_cohort = pd.read_csv(info_cohorts, sep="\t")
    # Add information of cancer type
    df_final_signatures_info = pd.merge(df_final_signatures, df_info_cohort)
    # Filter by expression
    df_final_signatures_info = filter_by_expression(df_final_signatures_info, exp_tcga)
    # Filter by Polymorphism
    df_final_signatures_info = filter_by_polymorphism(df_final_signatures_info, germline)
     # Add filter by olfactory receptors
    df_final_signatures_info = filter_by_olfactory_receptors(df_final_signatures_info, olfatory_receptors)
    # Add filter by known artifacts
    artifacts_file = path.join(FOLDER, 'data', "artifacts.json")
    with open(artifacts_file) as f:
        artifacts = json.load(f)
    df_final_signatures_info["Warning_Artifact"] = df_final_signatures_info.apply(lambda row: row["GENE"] in artifacts["suspects"], axis=1)
    df_final_signatures_info["Known_Artifact"] = df_final_signatures_info.apply(lambda row: row["GENE"] in artifacts["known"], axis=1)
    # Add filter by literature
    df_final_signatures_info = include_literature(df_final_signatures_info, cancermine)
    return df_final_signatures_info


def two_methods(vals):
    output = []
    for val in vals:
        if val is np.nan or str(val) == "nan":
            output.append(False)
            continue
        # Its a string
        list_values = val.split(",")
        if len(list_values) < 2:  # Only one value
            output.append(False)
            continue
        # Two methods
        methods = set(list_values)
        if len(methods) >= 2:
            output.append(True)
            continue
        output.append(False)
    return output


def get_cancer_genes(row):

    if row["CGC_GENE"] and (not row["CANCER_TYPE"] is np.nan) and str(row["CANCER_TYPE"]) in str(row["cancer_type_intogen"]):
        return True
    else:
        return False

def perform_vetting(df):
    l_data = []
    germ_center = ["AML","LY","CLL","MDS","DLBCL","NHLY"]
    for index, row in df.iterrows():
        if row["Warning_Expression"]:
            l = list(row.values)
            l.append("Warning expression")
            l_data.append(l)
        elif ((row["Signature9"] >= 0.5) and (row["CANCER_TYPE"] in germ_center)):
            l = list(row.values)
            l.append("Warning Signature9")
            l_data.append(l)
        elif row["Samples_3muts"] >= 1 and not(row["CGC_GENE"]):
            l = list(row.values)
            l.append("Samples with more than 3 mutations")
            l_data.append(l)
        elif row["MUTS/SAMPLE"] > 1.0 and row["Warning_Germline"] and not(row["Tier_CGC"]==1):
            l = list(row.values)
            l.append("Germline Warning")
            l_data.append(l)
        elif row["OR_Warning"]:
            l = list(row.values)
            l.append("Olfactory Receptor")
            l_data.append(l)
        elif row["Known_Artifact"]:
            l = list(row.values)
            l.append("Known artifact")
            l_data.append(l)
        elif row["n_papers"]== 0 and not(row["Tier_CGC"]==1):
            l = list(row.values)
            l.append("Lack of literature evidence")
            l_data.append(l)
        else:
            l = list(row.values)
            l.append("PASS")
            l_data.append(l)

    columns = list(df.columns) + ["FILTER"]
    df_filtered = pd.DataFrame(l_data, columns=columns)
    return df_filtered


def get_drivers(row):
    if row["TIER"] <= 3 and row["CGC_GENE"]:     # Tier 1 and tier 2 if cgc no bidders needed
        return True
    elif (len(str(row["Significant_Bidders"]).split(",")) > 1) and row["TIER"] <= 3 :     # Tier 1 if not cgc one bidder, #len(str(row["Significant_Bidders"]).split(",")) > 1str(row["Significant_Bidders"]) != "nan"
        return True
    else:
        return False


def concat(grp):
    set_ttypes = set(grp)
    return ",".join(list(set_ttypes))


def vet(df_vetting, paths, info_cohorts, cgc_path, threshold="05"):
    """Compute the driver list from the output of intogen and the vetting information"""

    l_data = []
    for path_ in paths:
        for file_data in glob.glob(path.join(path_, "*." + threshold + ".out.gz")):
            cohort = file_data.split("/")[-1].split(".")[0]
            df_data = pd.read_csv(file_data, sep="\t")
            df_data["COHORT"] = cohort
            df_data["PATH"] = file_data
            if not "PANEL" in cohort:  # Panels are not considered for the discovery
                l_data.append(df_data)
    # load cgc
    cgc = pd.read_csv(cgc_path, sep="\t")
    cgc["CGC_GENE"] = True
    cgc.rename(columns={"cancer_type": "cancer_type_intogen", "Tier": "Tier_CGC"}, inplace=True)
    # Read the data
    df_final = pd.concat(l_data, sort=True)
    # Load information of the data
    df_info = pd.read_csv(info_cohorts, sep="\t")
    df_info.rename(columns={"MUTATIONS": "MUTATIONS_COHORT", "SAMPLES": "SAMPLES_COHORT"}, inplace=True)
    df = df_final.merge(df_info, how="left", left_on="COHORT", right_on="COHORT")
    df = pd.merge(df, cgc[["Gene Symbol", "CGC_GENE", "cancer_type_intogen", "Tier_CGC"]], left_on="SYMBOL",
                  right_on="Gene Symbol", how="left")
    df["CGC_GENE"].fillna(False, inplace=True)
    df["driver"] = df.apply(lambda row: get_drivers(row), axis=1)
    df_drivers = df[df["driver"]]
    print("Number of drivers pre-vetting:" + str(len(df_drivers["SYMBOL"].unique())))
    # Include cgc
    df_drivers["CGC_CANCER_GENE"] = df_drivers.apply(lambda row: get_cancer_genes(row), axis=1)
    df_drivers.drop(["Gene Symbol", "cancer_type_intogen"], inplace=True, axis=1)
    # Include average number of mutations per sample
    df_drivers["MUTS/SAMPLE"] = df_drivers.apply(lambda row: row["MUTS"] / row["SAMPLES"], axis=1)
    # Include the number of cohorts per gene
    # Add warning of number of cohorts per gene
    df_counts = df_drivers.groupby("SYMBOL", as_index=False).agg({"COHORT": "count"})
    df_counts.rename(columns={"COHORT": "num_cohorts"}, inplace=True)
    df_drivers = df_drivers.merge(df_counts)
    df_drivers["Warning_num_cohorts"] = df_drivers.apply(lambda row: True if row["num_cohorts"] == 1 else False, axis=1)

    # Perform the vetting
    df_vetting.rename(columns={"GENE": "SYMBOL"}, inplace=True)
    df_drivers_vetting = pd.merge(df_drivers, df_vetting[
        ["SNP", "INDEL", "COHORT", "INDEL/SNP", "Signature10", "Signature9", "Warning_Expression", "Warning_Germline",
         "SYMBOL", "Samples_3muts", "OR_Warning", "Warning_Artifact", "Known_Artifact", "n_papers"]].drop_duplicates(), how="left")
    df_drivers_vetting["Warning_Expression"].fillna(False, inplace=True)
    df_drivers_vetting["Warning_Germline"].fillna(False, inplace=True)
    df_drivers_vetting["OR_Warning"].fillna(False, inplace=True)
    df_drivers_vetting["Warning_Artifact"].fillna(False, inplace=True)
    df_drivers_vetting["Known_Artifact"].fillna(False, inplace=True)
    df_drivers_vetting["Signature9"].fillna(0.0, inplace=True)
    df_drivers_vetting["Signature10"].fillna(0.0, inplace=True)
    df_drivers_vetting["Samples_3muts"].fillna(0.0, inplace=True)
    df_drivers_vetting_info = perform_vetting(df_drivers_vetting)
    print("Number of drivers after-vetting:" + str(
        len(df_drivers_vetting_info[df_drivers_vetting_info["FILTER"] == "PASS"]["SYMBOL"].unique())))
    return df_drivers_vetting_info


def read_file(filein):
    f = open(filein, 'r')
    genes = set()
    for line in f.readlines():
        line = line.strip()
        genes.add(line)
    f.close()
    return genes


def filter(df_all):

    # Read vetted file
    df_drivers = df_all[df_all["FILTER"] == "PASS"]
    print("Number of drivers before white/black listing:" + str(len(df_drivers["SYMBOL"].unique())))

    # Remove black listed genes
    black_listed_file = path.join(FOLDER, 'data', 'black_listed.txt')
    black_listed = read_file(black_listed_file)
    df_drivers = df_drivers[~df_drivers["SYMBOL"].isin(black_listed)]

    # Now rescue white listed discarded genes
    df_discarded = df_all[df_all["FILTER"] == "Lack of literature evidence"]

    # Recover white listed genes
    white_listed_file = path.join(FOLDER, 'data', 'white_listed.txt')
    white_listed = read_file(white_listed_file)
    df_recovered = df_discarded[df_discarded["SYMBOL"].isin(white_listed)]
    df_recovered["FILTER"] = "PASS"

    df_final_list = pd.concat([df_drivers, df_recovered], sort=True)
    print("Number of drivers after white/black listing:" + str(
        len(df_final_list[df_final_list["FILTER"] == "PASS"]["SYMBOL"].unique())))

    return df_final_list


def run(paths, expression, olfactory, cancermine, germline, cgc, ensembl, cohorts,
        output, tmp_folder=None, muts=3, threshold="05"):
    tmp_folder = tmp_folder or tempfile.TemporaryDirectory().name

    preparation_paths = [path.join(path_, 'deconstructsig') for path_ in paths]
    df = prepare(preparation_paths, expression, olfactory, cancermine, germline, cohorts, muts)
    prepared_file = path.join(tmp_folder, 'information_vetting_genes.tsv.gz')
    df.to_csv(prepared_file, sep="\t", index=False)

    vet_paths = [path.join(path_, 'combination') for path_ in paths]
    df = vet(df, vet_paths, cohorts, cgc, threshold)
    vet_file = path.join(tmp_folder, f'all_drivers{threshold}.tsv')
    df.rename(columns={"QVALUE_stouffer_w":"QVALUE_COMBINATION","Significant_Bidders":"METHODS","n_papers":"NUM_PAPERS", "cancer_type":"CANCER_TYPE"},inplace=True)
    df.columns = map(str.upper, df.columns)
    columns= ["SYMBOL","COHORT","CANCER_TYPE","METHODS","QVALUE_COMBINATION","TIER","MUTATIONS_COHORT","SAMPLES_COHORT","ROLE", "CGC_GENE", "TIER_CGC",
              "CGC_CANCER_GENE", "NUM_COHORTS", "SIGNATURE9", "WARNING_EXPRESSION", "WARNING_GERMLINE", "SAMPLES_3MUTS", "OR_WARNING",
              "KNOWN_ARTIFACT", "NUM_PAPERS", "FILTER"]
    df[columns].sort_values(["SYMBOL","CANCER_TYPE"]).to_csv(vet_file, sep="\t", index=False)

    df = filter(df)

    df_biomart = pd.read_csv(ensembl, sep="\t", index_col=False, usecols=[1, 10],
                             names=["SYMBOL", "TRANSCRIPT"],
                             header=None)
    df_biomart.drop_duplicates(inplace=True)
    duplicated = df_biomart[df_biomart.duplicated(subset='SYMBOL', keep=False)]['SYMBOL'].unique()

    symbols = df['SYMBOL'].values
    if any(x in symbols for x in duplicated):
        logging.error('A CGC symbol appears mapped to 2+ transcripts')
        quit(-1)

    df.to_csv(output, sep='\t', index=False)

    # Create a unique file of drivers
    drivers=df[df["FILTER"] == "PASS"]["SYMBOL"].unique()
    # Add the ensembl gene id
    df_ensembl = pd.read_csv(ensembl, sep="\t", index_col=False, usecols=[0,1,2,10], names=["ENSEMBL_GENE","SYMBOL","ENSEMBL_PROTEIN","ENSEMBL_TRANSCRIPT"], header=None) # ENSG00000160752	FDPS	ENSP00000349078	1	155312255	155312395	340	480	1260	1	ENST00000356657
    df_drivers_unique= df_ensembl[df_ensembl["SYMBOL"].isin(drivers)].drop_duplicates()
    unique_drivers_file = path.join(tmp_folder, f'unique_drivers{threshold}.tsv')
    df_drivers_unique.to_csv(unique_drivers_file, sep="\t",index=False)


if __name__ == "__main__":
    output = sys.argv[1]
    expression = sys.argv[2]
    orfactory = sys.argv[3]
    cancermine = sys.argv[4]
    germline = sys.argv[5]
    cgc = sys.argv[6]
    ensembl = sys.argv[7]
    cohorts = sys.argv[8]
    threshold = sys.argv[9]
    paths = sys.argv[10:-1]
    tmp_folder = sys.argv[-1]
    run(paths, expression, orfactory, cancermine, germline, cgc, ensembl, cohorts,
        output, tmp_folder, threshold=threshold)
