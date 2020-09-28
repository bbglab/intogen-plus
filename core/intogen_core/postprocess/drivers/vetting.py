
import os

import pandas as pd


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


def perform_vetting(df, ctype):
    l_data = []
    germ_center = ["AML","LY","CLL","MDS","DLBCL","NHLY"]
    for index, row in df.iterrows():
        if row["Warning_Expression"]:
            l = list(row.values)
            l.append("Warning expression")
            l_data.append(l)
        elif ((row["Signature9"] >= 0.5) and (ctype in germ_center)):
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


def vet(df_vetting, combination, ctype):
    """Compute the driver list from the output of intogen and the vetting information"""

    df = pd.read_csv(combination, sep="\t")

    # load cgc
    cgc_path = os.path.join(os.environ['INTOGEN_DATASETS'], 'cgc', 'cancer_gene_census_parsed.tsv')
    cgc = pd.read_csv(cgc_path, sep="\t")
    cgc["CGC_GENE"] = True
    cgc.rename(columns={"cancer_type": "cancer_type_intogen", "Tier": "Tier_CGC"}, inplace=True)

    df = pd.merge(df, cgc[["Gene Symbol", "CGC_GENE", "cancer_type_intogen", "Tier_CGC"]],
                  left_on="SYMBOL", right_on="Gene Symbol", how="left")
    df["CGC_GENE"].fillna(False, inplace=True)
    if len(df) == 0:
        # Simply add the column
        df['driver'] = None
    else:
        df["driver"] = df.apply(lambda row: get_drivers(row), axis=1)
    df_drivers = df[df["driver"]]
    print("Number of drivers pre-vetting:" + str(len(df_drivers["SYMBOL"].unique())))

    if len(df_drivers["SYMBOL"].unique()) == 0:
        # Simply add the columns
        df_drivers["CGC_CANCER_GENE"] = None
        df_drivers["MUTS/SAMPLE"] = None
    else:
        # Include cgc
        df_drivers["CGC_CANCER_GENE"] = df_drivers.apply(lambda row: get_cancer_genes(row, ctype), axis=1)
        df_drivers.drop(["Gene Symbol", "cancer_type_intogen"], inplace=True, axis=1)
        # Include average number of mutations per sample
        df_drivers["MUTS/SAMPLE"] = df_drivers.apply(lambda row: row["MUTS"] / row["SAMPLES"], axis=1)
        # Include the number of cohorts per gene

    # Perform the vetting
    df_vetting.rename(columns={"GENE": "SYMBOL"}, inplace=True)
    df_drivers_vetting = pd.merge(df_drivers, df_vetting[
        ["SNP", "INDEL", "INDEL/SNP", "Signature10", "Signature9", "Warning_Expression", "Warning_Germline",
         "SYMBOL", "Samples_3muts", "OR_Warning", "Warning_Artifact", "Known_Artifact", "n_papers"]].drop_duplicates(), how="left")
    df_drivers_vetting["Warning_Expression"].fillna(False, inplace=True)
    df_drivers_vetting["Warning_Germline"].fillna(False, inplace=True)
    df_drivers_vetting["OR_Warning"].fillna(False, inplace=True)
    df_drivers_vetting["Warning_Artifact"].fillna(False, inplace=True)
    df_drivers_vetting["Known_Artifact"].fillna(False, inplace=True)
    df_drivers_vetting["Signature9"].fillna(0.0, inplace=True)
    df_drivers_vetting["Signature10"].fillna(0.0, inplace=True)
    df_drivers_vetting["Samples_3muts"].fillna(0.0, inplace=True)
    df_drivers_vetting_info = perform_vetting(df_drivers_vetting, ctype)
    print("Number of drivers after-vetting:" + str(
        len(df_drivers_vetting_info[df_drivers_vetting_info["FILTER"] == "PASS"]["SYMBOL"].unique())))

    return df_drivers_vetting_info
