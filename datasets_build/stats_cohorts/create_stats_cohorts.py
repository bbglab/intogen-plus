import pandas as pd
import os
import click
import json
import bglogs



def select_ttype(row,d_names):
    '''
    Return the ttype of the cohort
    :param row:
    :return:
    '''
    if "ICGC" in row["source"]:
        return d_names["_".join(row["COHORT"].split("_")[2:])]
    elif "PCAWG" in row["source"]:
        return d_names["_".join(row["COHORT"].split("_")[2:])]
    elif "HARTWIG" in row["source"]:
        return d_names["_".join(row["COHORT"].split("_")[2:])]
    elif "TCGA" in row["source"]:
        return d_names[row["COHORT"].split("_")[2]]
    elif "CBIOP" in row["source"]:
        return d_names[row["COHORT"].split("_")[2]]
    elif "OTHER" in row["source"]:
        return d_names["_".join(row["COHORT"].split("_")[2:])]
    elif "PEDCBIOP" in row["source"]:
        return d_names[row["COHORT"].split("_")[2]]
    elif "ST_JUDE" in row["source"]:
        return d_names[row["COHORT"].split("_")[2]]
    elif "TARGET" in row["source"]:
        return d_names[row["COHORT"].split("_")[2]]
    else:
        return np.nan

def read_stjude(path_stjude):
    '''
    Read the datasets from stjude and return a Dataframe with all the information of the cohorts
    :param path_stjude:
    :return:
    '''
    df_st_jude = pd.read_csv(os.path.join(path_stjude, "cancertypes_mapping.txt"), sep="\t",
                             names=["ORIGINAL_NAME", "NAME", "cancer_type"])
    df_samples_stjude = pd.read_csv(os.path.join(path_stjude, "count_samples.txt"), names=["SAMPLES", "NAME"],
                                    sep=" ")
    df_vars_stjude = pd.read_csv(os.path.join(path_stjude, "count_variants.txt"), names=["NAME", "MUTATIONS"],
                                 sep="\t")
    df_data_stjude = df_samples_stjude.merge(df_vars_stjude)
    df_data_stjude = df_data_stjude.merge(df_st_jude[["NAME", "cancer_type"]].drop_duplicates())

    dict_ttypes_stjude = {"M": "Metastasis", "D": "Primary", "R": "Relapse"}
    df_data_stjude["source"] = "STJUDE"
    df_data_stjude["PLATFORM"] = "WGS"
    df_data_stjude["COHORT"] = df_data_stjude["NAME"]
    df_data_stjude["TYPE"] = df_data_stjude.apply(lambda row: dict_ttypes_stjude[row["COHORT"].split("_")[0]], axis=1)
    df_data_stjude["AGE"] = "Pediatric"
    df_data_stjude["TREATED"] = df_data_stjude.apply(lambda row: "Treated" if "R_" in row["COHORT"] or "M_" in row["COHORT"] else "Untreated", axis=1)
    df_data_stjude.drop(columns=["NAME"], inplace=True)
    return df_data_stjude

def read_hartwig(path_hartwig,d_names):
    '''
    Read the datasets from hartwig and return a Dataframe with all the information of the cohorts
    :param path_hartwig:
    :return:
    '''
    df_samples_hartwig = pd.read_csv(os.path.join(path_hartwig, "count_samples.txt"), names=["SAMPLES", "NAME"], sep=" ")
    df_vars_hartiwg = pd.read_csv(os.path.join(path_hartwig, "count_variants.txt"), names=["NAME", "MUTATIONS"], sep="\t")
    df_data_hartwig = df_samples_hartwig.merge(df_vars_hartiwg)
    df_data_hartwig["source"] = "HARTWIG"
    df_data_hartwig["PLATFORM"] = "WGS"
    df_data_hartwig["cancer_type"] = df_data_hartwig.apply(
        lambda row: d_names[row["NAME"]] if row["NAME"] in d_names else row["NAME"], axis=1)
    df_data_hartwig["COHORT"] = df_data_hartwig["NAME"]
    df_data_hartwig["TYPE"] = "Metastasis"
    df_data_hartwig["AGE"] = "Adult"
    df_data_hartwig["TREATED"] = "Treated"
    df_data_hartwig.drop(columns=["NAME"], inplace=True)
    return df_data_hartwig

def read_intogen(path_intogen,d_names):
    '''
    Read the datasets from intogen and return a Dataframe with all the information of the cohorts
    :param path_intogen:
    :param d_names:
    :return:
    '''
    df_samples_intogen = pd.read_csv(os.path.join(path_intogen, "count_samples.txt"), names=["SAMPLES", "COHORT"],
                                  sep=" ")

    df_vars_intogen = pd.read_csv(os.path.join(path_intogen, "count_variants.txt"), names=["COHORT", "MUTATIONS"], sep="\t")
    df_data_intogen = df_samples_intogen.merge(df_vars_intogen)
    df_data_intogen["source"] = df_data_intogen.apply(lambda row: row['COHORT'].split('_')[0], axis=1)
    path_info_extra = os.path.join(path_intogen, "info_datasets.csv")
    df_info_extra = pd.read_csv(path_info_extra, sep="\t")
    df_data_intogen = df_data_intogen.merge(
    df_info_extra[["COHORT", "PLATFORM", "TYPE", "AGE", "TREATED"]].drop_duplicates(), how="left")
    df_data_intogen["cancer_type"] = df_data_intogen.apply(lambda row: select_ttype(row,d_names), axis=1)
    return df_data_intogen


def read_data(path_harwig,path_intogen,path_stjude,dict_ttypes):
    # Read dictionary of ttypes per cohort
    with open(dict_ttypes, 'r') as fp:
        d_names = json.load(fp)
    # Read st jude
    df_stjude = read_stjude(path_stjude)
    # Read hartwig
    df_harwig = read_hartwig(path_harwig,d_names)
    # Read intogen
    df_intogen = read_intogen(path_intogen,d_names)
    # Concat all information
    df_final = pd.concat([df_intogen, df_harwig, df_stjude], sort=True)

    # Include the name of the tumor type
    with open("dictionary_labels.json",'r') as f:
        d_legend = json.load(f)
    df_final["legend"] = df_final.apply(
        lambda row: d_legend[row["cancer_type"]] if row["cancer_type"] in d_legend else row["cancer_type"],
        axis=1)
    # Include the long name of tumor type
    with open("dictionary_long_names.json",'r') as f:
        d_name_web = json.load(f)
    df_final["web_name"] = df_final.apply(
        lambda row: d_name_web[row["cancer_type"]] if row["cancer_type"] in d_name_web else row["cancer_type"],
        axis=1)
    return df_final



@click.command()
@click.option('--path_hartwig',help= 'path to the datasets from Hartwig', type=click.Path(),required=True)
@click.option('--path_intogen',help= 'path to the datasets from intOGen', type=click.Path(),required=True)
@click.option('--path_stjude',help= 'path to the datasets from stjude', type=click.Path(),required=True)
@click.option('--dict_ttypes',help= 'mapping cohorts to tumor type', type=click.Path(),required=True)
@click.option('--path_output',help= 'output file to store the results', type=click.Path(),required=True)
@click.option('--debug', is_flag=True)
def cmdline(path_hartwig, path_intogen, path_stjude, dict_ttypes, path_output, debug):
    bglogs.configure(debug=debug)
    df = read_data(path_hartwig,path_intogen,path_stjude,dict_ttypes)
    bglogs.info("Storing results in {}",path_output)
    df.to_csv(path_output,sep="\t",index=False)

if __name__ == "__main__":
    cmdline()