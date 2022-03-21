
import os

import pandas as pd

from intogen_core.exceptions import IntogenError
from intogen_core.postprocess.drivers.bw_list import read_file


def get_drivers(row):
    if row["TIER"] <= 3 and row["CGC_GENE"]:     # Tier 1 and tier 2 if cgc no bidders needed
        return True
    elif (len(str(row["Significant_Bidders"]).split(",")) > 1) and row["TIER"] <= 3 :     # Tier 1 if not cgc one bidder, #len(str(row["Significant_Bidders"]).split(",")) > 1str(row["Significant_Bidders"]) != "nan"
        return True
    else:
        return False


def get_cancer_genes(row, ctype):

    if row["CGC_GENE"] and ctype is not None and ctype in str(row["cancer_type_intogen"]):
        return True
    else:
        return False


def perform_vetting(row, ctype):

    germ_center = ["LYMPHOID","CLLSLL","DLBCLNOS","NHL","ALL","LNM","TLL","BLL","HL"] #all lymphoid cancer types

    #Add white list for literature vetting; 
    white_listed_file = os.path.join(os.environ['INTOGEN_DATASETS'], 'postprocess', 'white_listed.txt')
    white_listed = read_file(white_listed_file)
    
    if row['DRIVER'] == False:
        out = "No driver"
    elif row["Warning_Expression"]:
        out ="Warning expression"
    elif ((row["Signature9"] >= 0.5) and (ctype in germ_center)):
        out = "Warning Signature9"
    elif row["Samples_3muts"] >= 1 and not(row["CGC_GENE"]):
        out = "Samples with more than 2 mutations"
    elif row["MUTS/SAMPLE"] > 1.0 and row["Warning_Germline"] and not(row["CGC_GENE"]):
        out = "Germline Warning"
    elif row["OR_Warning"]:
        out = "Olfactory Receptor"
    elif row["Known_Artifact"]:
        out = "Known artifact"
    elif row["Warning_Artifact"]:
        out = "Warning artifact"
    elif row["n_papers"]== 0 and not(row["CGC_GENE"]) and (row['SYMBOL'] not in white_listed):
        out = "Lack of literature evidence"
    else:
        out = "PASS"
    return out

def vet(df_vetting, combination, ctype):
    """Compute the driver list from the output of intogen and the vetting information"""

    df = pd.read_csv(combination, sep="\t")
    if len(df) == 0:
        #raise IntogenError('No drivers in combination to perform vetting')
        print('No drivers in combination to perform vetting')
        df_drivers_vetting = pd.DataFrame()
    else:
        # load cgc
        cgc_path = os.path.join(os.environ['INTOGEN_DATASETS'], 'cgc', 'cancer_gene_census_parsed.tsv')
        cgc = pd.read_csv(cgc_path, sep="\t")
        cgc["CGC_GENE"] = True
        cgc.rename(columns={"cancer_type": "cancer_type_intogen", "Tier": "Tier_CGC"}, inplace=True)

        df_drivers = pd.merge(df, cgc[["Gene Symbol", "CGC_GENE", "cancer_type_intogen", "Tier_CGC"]],
                      left_on="SYMBOL", right_on="Gene Symbol", how="left")
        df_drivers["CGC_GENE"].fillna(False, inplace=True)
        df_drivers["DRIVER"] = df_drivers.apply(lambda row: get_drivers(row), axis=1)
        print("Number of drivers pre-vetting:" + str(len(df_drivers["SYMBOL"][df_drivers["DRIVER"]==True].unique())))

        if len(df_drivers["SYMBOL"][df_drivers["DRIVER"]==True].unique()) == 0:
            # Simply add the columns
            df_drivers["CGC_CANCER_GENE"] = None
            df_drivers["MUTS/SAMPLE"] = None
        else:
            # Include cgc
            df_drivers["CGC_CANCER_GENE"] = df_drivers.apply(lambda row: get_cancer_genes(row, ctype), axis=1)
            df_drivers.drop(["Gene Symbol", "cancer_type_intogen"], inplace=True, axis=1)
            # Include average number of mutations per sample
            df_drivers["MUTS/SAMPLE"] = df_drivers.apply(lambda row: row["MUTS"] / row["SAMPLES"], axis=1)
            # Include the number of cohorts per gene (?)

        # Perform the vetting
        df_vetting.rename(columns={"GENE": "SYMBOL"}, inplace=True)
        df_drivers_vetting = pd.merge(df_drivers, df_vetting[
            ["SNV", "INDEL", "INDEL/SNV", "Signature10", "Signature9", "Warning_Expression", "Warning_Germline",
             "SYMBOL", "Samples_3muts", "OR_Warning", "Warning_Artifact", "Known_Artifact", "n_papers"]].drop_duplicates(), how="left")
        df_drivers_vetting["DRIVER"].fillna(False, inplace=True)
        df_drivers_vetting["Warning_Expression"].fillna(False, inplace=True)
        df_drivers_vetting["Warning_Germline"].fillna(False, inplace=True)
        df_drivers_vetting["OR_Warning"].fillna(False, inplace=True)
        df_drivers_vetting["Warning_Artifact"].fillna(False, inplace=True)
        df_drivers_vetting["Known_Artifact"].fillna(False, inplace=True)
        df_drivers_vetting["Signature9"].fillna(0.0, inplace=True)
        df_drivers_vetting["Signature10"].fillna(0.0, inplace=True)
        df_drivers_vetting["Samples_3muts"].fillna(0.0, inplace=True)
        df_drivers_vetting["FILTER"] = df_drivers_vetting.apply(lambda row: perform_vetting(row, ctype), axis=1)
        print("Number of drivers after-vetting:" + str(
            len(df_drivers_vetting[df_drivers_vetting["FILTER"] == "PASS"]["SYMBOL"].unique())))

    return df_drivers_vetting
