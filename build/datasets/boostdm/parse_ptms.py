import click
import pandas as pd
import os
import json

def get_aa(row):
    return row["MOD_RSD"].split("-")[0]

def read_data(base_path):
    # Acetylation
    acetilation = pd.read_csv(os.path.join(base_path, "Acetylation_site_dataset.gz"), sep="\t", compression="gzip",
                              skiprows=2)
    acetilation = acetilation[acetilation["ORGANISM"] == "human"]
    acetilation["Acetylation"] = 1
    acetilation["AA"] = acetilation.apply(lambda row: get_aa(row), axis=1)
    acetilation = acetilation[["GENE", "AA", "Acetylation"]].drop_duplicates()
    # Glycosylation
    galnac = pd.read_csv(os.path.join(base_path, "O-GalNAc_site_dataset.gz"), sep="\t", compression="gzip", skiprows=2)
    galnac = galnac[galnac["ORGANISM"] == "human"]
    galnac["O-glycosylation"] = 1
    galnac["AA"] = galnac.apply(lambda row: get_aa(row), axis=1)
    galnac = galnac[["GENE", "AA", "O-glycosylation"]].drop_duplicates()
    # GlcNAc
    glcnac = pd.read_csv(os.path.join(base_path, "O-GlcNAc_site_dataset.gz"), sep="\t", compression="gzip", skiprows=2)
    glcnac = glcnac[glcnac["ORGANISM"] == "human"]
    glcnac["O-GlcNAc"] = 1
    glcnac["AA"] = glcnac.apply(lambda row: get_aa(row), axis=1)
    glcnac = glcnac[["GENE", "AA", "O-GlcNAc"]].drop_duplicates()
    # Phosphorylation
    phospho = pd.read_csv(os.path.join(base_path, "Phosphorylation_site_dataset.gz"), sep="\t", compression="gzip",
                          skiprows=2)
    phospho = phospho[phospho["ORGANISM"] == "human"]
    phospho["Phosphorylation"] = 1
    phospho["AA"] = phospho.apply(lambda row: get_aa(row), axis=1)
    phospho = phospho[["GENE", "AA", "Phosphorylation"]].drop_duplicates()
    # Ubiquitination
    ubi = pd.read_csv(os.path.join(base_path, "Ubiquitination_site_dataset.gz"), sep="\t", compression="gzip",
                      skiprows=2)
    ubi = ubi[ubi["ORGANISM"] == "human"]
    ubi["Ubiquitination"] = 1
    ubi["AA"] = ubi.apply(lambda row: get_aa(row), axis=1)
    ubi = ubi[["GENE", "AA", "Ubiquitination"]].drop_duplicates()
    # Methylation
    methyl = pd.read_csv(os.path.join(base_path, "Methylation_site_dataset.gz"), sep="\t", compression="gzip",
                         skiprows=2)
    methyl = methyl[methyl["ORGANISM"] == "human"]
    methyl["Methylation"] = 1
    methyl["AA"] = methyl.apply(lambda row: get_aa(row), axis=1)
    methyl = methyl[["GENE", "AA", "Methylation"]].drop_duplicates()
    # Sumoylation
    sumo = pd.read_csv(os.path.join(base_path, "Sumoylation_site_dataset.gz"), sep="\t", compression="gzip", skiprows=2)
    sumo = sumo[sumo["ORGANISM"] == "human"]
    sumo["Sumoylation"] = 1
    sumo["AA"] = sumo.apply(lambda row: get_aa(row), axis=1)
    sumo = sumo[["GENE", "AA", "Sumoylation"]].drop_duplicates()
    # Regulatory sites
    regulatory = pd.read_csv(os.path.join(base_path, "Regulatory_sites.gz"), sep="\t", compression="gzip", skiprows=3,
                             error_bad_lines=False)
    regulatory = regulatory[regulatory["ORGANISM"] == "human"]
    regulatory["Regulatory_Site"] = 1
    regulatory["AA"] = regulatory.apply(lambda row: get_aa(row), axis=1)
    regulatory = regulatory[["GENE", "AA", "Regulatory_Site"]].drop_duplicates()
    list_data = [acetilation, galnac, glcnac, phospho, ubi, methyl, sumo, regulatory]
    # Create dataframe with information
    df_ptms = list_data[0]
    for d in list_data[1:]:
        df_ptms = df_ptms.merge(d, how="outer")
    # Save it as a dictionary with key gene and amino acid
    df_ptms.fillna(0, inplace=True)
    d_data = {}
    for index, row in df_ptms.iterrows():
        if not (row["GENE"] in d_data):
            d_data[row["GENE"]] = {}
        d_data[row["GENE"]][row["AA"]] = [row['Acetylation'], row['O-glycosylation'], row['O-GlcNAc'],
                                          row['Phosphorylation'], row['Ubiquitination'], row['Methylation'],
                                          row['Sumoylation'], row['Regulatory_Site']]
    return d_data


@click.command()
@click.option('--path_raw_data',help= 'path to the raw datasets from PhosphositePlus', type=click.Path(),required=True)
@click.option('--path_output_dictionary',help= 'path to the output .json file', type=click.Path(),required=True)
def cmdline(path_raw_data, path_output_dictionary):
    dictionary = read_data(path_raw_data)
    with open(path_output_dictionary,"w") as write_file:
        json.dump(dictionary, write_file)

if __name__ == "__main__":
    cmdline()


