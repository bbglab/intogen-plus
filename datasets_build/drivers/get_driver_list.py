'''
Compute the driver list from the output of intogen and the vetting information
'''

import pandas as pd
pd.options.mode.chained_assignment = None
import glob
import json
import os
import numpy as np
import json
import click



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

    if row["CGC_GENE"] and (not row["cancer_type"] is np.nan) and str(row["cancer_type"]) in str(row["cancer_type_intogen"]):
        return True
    else:
        return False

def perform_vetting(df):
    with open("artifacts.json") as f:
        artifacts = json.load(f)
    l_data = []
    germ_center = ["AML","LY","CLL","MDS","DLBCL","NHLY"]
    for index, row in df.iterrows():
        if row["Warning_Expression"]:
            l = list(row.values)
            l.append("Warning expression")
            l_data.append(l)
        elif ((row["Signature9"] >= 0.5) and (row["cancer_type"] in germ_center)):
            l = list(row.values)
            l.append("Warning Signature9")
            l_data.append(l)
        elif row["Warning_num_cohorts"] and not(row["CGC_GENE"]) and (row["Warning_Artifact"]):
            l = list(row.values)
            l.append("Warning single cohort and lack evidence")
            l_data.append(l)
        elif row["Samples_3muts"] >= 1 and not(row["CGC_GENE"]):
            l = list(row.values)
            l.append("Samples with more than 3 mutations")
            l_data.append(l)
        elif row["MUTS/SAMPLE"] > 1.0 and row["Warning_Germline"]: # Less than 5 samples mutated and warning of germline
            l = list(row.values)
            l.append("Germline Warning")
            l_data.append(l)
        elif row["OR_Warning"]:
            l = list(row.values)
            l.append("Olfactory Receptor")
            l_data.append(l)
        elif row["Warning_Artifact"] or row["SYMBOL"] in artifacts["suspects"]:
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


def main(paths,info_cohorts,dir_out,threshold,cgc_path,vetting_file,ensembl_file):
    l_data = []
    for path in paths:
        for file_data in glob.glob(os.path.join(path, "*." + threshold + ".out.gz")):
            cohort = file_data.split("/")[-1].split(".")[0]
            df_data = pd.read_csv(file_data, sep="\t")
            df_data["COHORT"] = cohort
            df_data["PATH"] = file_data
            if not "PANEL" in cohort: # Panels are not considered for the discovery
                l_data.append(df_data)
    # load cgc
    cgc = pd.read_csv(cgc_path, sep="\t")
    cgc["CGC_GENE"] = True
    cgc.rename(columns={"cancer_type":"cancer_type_intogen"},inplace=True)
    # Read the data
    df_final = pd.concat(l_data,sort=True)
    # Load information of the data
    df_info = pd.read_csv(info_cohorts,sep="\t")
    df_info.rename(columns={"MUTATIONS": "MUTATIONS_COHORT", "SAMPLES": "SAMPLES_COHORT"}, inplace=True)
    df = df_final.merge(df_info, how="left", left_on="COHORT", right_on="COHORT")
    df = pd.merge(df, cgc[["Gene Symbol", "CGC_GENE", "cancer_type_intogen"]], left_on="SYMBOL",
             right_on="Gene Symbol", how="left")
    df["CGC_GENE"].fillna(False, inplace=True)
    df["driver"] = df.apply(lambda row: get_drivers(row),axis=1)
    df_drivers = df[df["driver"]]
    print ("Number of drivers pre-vetting:" + str(len(df_drivers["SYMBOL"].unique())))
    # Include cgc
    df_drivers["CGC_CANCER_GENE"] = df_drivers.apply(lambda row: get_cancer_genes(row), axis=1)
    df_drivers.drop(["Gene Symbol", "cancer_type_intogen"], inplace=True, axis=1)
    # Include average number of mutations per sample
    df_drivers["MUTS/SAMPLE"] = df_drivers.apply(lambda row: row["MUTS"] / row["SAMPLES"], axis=1)
    # Include the number of cohorts per gene
    # Add warning of number of cohorts per gene
    df_counts = df_drivers.groupby("SYMBOL",as_index=False).agg({"COHORT":"count"})
    df_counts.rename(columns={"COHORT":"num_cohorts"},inplace=True)
    df_drivers=df_drivers.merge(df_counts)
    df_drivers["Warning_num_cohorts"] = df_drivers.apply(lambda row: True if row["num_cohorts"] == 1 else False,axis=1)

    # Perform the vetting
    df_vetting = pd.read_csv(vetting_file, sep="\t",
                             compression="gzip")

    df_vetting.rename(columns={"GENE": "SYMBOL"}, inplace=True)
    df_drivers_vetting = pd.merge(df_drivers, df_vetting[
        ["SNP", "INDEL", "COHORT", "INDEL/SNP", "Signature10", "Signature9", "Warning_Expression", "Warning_Germline",
        "SYMBOL", "Samples_3muts","OR_Warning","Warning_Artifact"]].drop_duplicates(), how="left")
    df_drivers_vetting["Warning_Expression"].fillna(False, inplace=True)
    df_drivers_vetting["Warning_Germline"].fillna(False, inplace=True)
    df_drivers_vetting["OR_Warning"].fillna(False, inplace=True)
    df_drivers_vetting["Warning_Artifact"].fillna(False, inplace=True)
    df_drivers_vetting["Signature9"].fillna(0.0, inplace=True)
    df_drivers_vetting["Signature10"].fillna(0.0, inplace=True)
    df_drivers_vetting["Samples_3muts"].fillna(0.0, inplace=True)
    df_drivers_vetting_info = perform_vetting(df_drivers_vetting)
    # Save all information
    df_drivers_vetting_info.to_csv(
        os.path.join(dir_out,"all_drivers"+threshold+".tsv"), sep="\t",index=False)
    # Save only those non-vetted
    df_drivers_vetting_info[df_drivers_vetting_info["FILTER"] == "PASS"].to_csv(os.path.join(dir_out,"vetted_drivers"+threshold+".tsv"), sep="\t",index=False)
    print ("Number of drivers after-vetting:" + str(len(df_drivers_vetting_info[df_drivers_vetting_info["FILTER"]=="PASS"]["SYMBOL"].unique())))
    # Create a unique file of drivers
    drivers=df_drivers_vetting_info[df_drivers_vetting_info["FILTER"] == "PASS"]["SYMBOL"].unique()
    # Add the ensembl gene id
    df_ensembl = pd.read_csv(ensembl_file, sep="\t", index_col=False, usecols=[0,1,2,10], names=["ENSEMBL_GENE","SYMBOL","ENSEMBL_PROTEIN","ENSEMBL_TRANSCRIPT"], header=None) # ENSG00000160752	FDPS	ENSP00000349078	1	155312255	155312395	340	480	1260	1	ENST00000356657
    df_drivers_unique= df_ensembl[df_ensembl["SYMBOL"].isin(drivers)].drop_duplicates()
    df_drivers_unique.to_csv(os.path.join(dir_out,"unique_drivers"+threshold+".tsv"), sep="\t",index=False)


@click.command()
@click.option('-i', '--intogen', 'intogen',type=click.Path(),  help='Path to intogen data',required=True)
@click.option('-w', '--hartwig', 'hartwig', type=click.Path(),  help='Path to hartwig data', required=True)
@click.option('-s', '--stjude', 'stjude', type=click.Path(),  help="Path to stjude data", required=True)
@click.option('-c', '--info_cohorts', 'info_cohorts', type=click.Path(),  help="Path to info of the cohorts", required=True)
@click.option('-o', '--output', 'output', type=click.Path(),  help="Path to the output folder.", required=True)
@click.option('-t', '--threshold', 'threshold',  help="Threshold of the output files, 05 or 01. Default 05", default="05")
@click.option('-g', '--cgc_path', 'cgc_path', type=click.Path(),  help="Path to the CGC information data")
@click.option('-n', '--ensembl_file', 'ensembl_file', type=click.Path(),  help="Path to the canonical trancript information file",default="/workspace/projects/intogen_2017/pipeline/datasets/hg38_vep92_develop/shared/cds_biomart.tsv")
@click.option('-v', '--vetting', 'vetting_file', type=click.Path(),  help="Path to the vetting file")

def cmdline(intogen, hartwig, stjude, info_cohorts, output, threshold, cgc_path, vetting_file,ensembl_file):
    paths = [intogen,hartwig,stjude]
    main(paths,info_cohorts,output,threshold,cgc_path,vetting_file,ensembl_file)



if __name__ == "__main__":
    cmdline()

