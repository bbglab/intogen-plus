import math
import os

import bgdata
import numpy as np
import pandas as pd


def excess_muts(n_obs, omega):
    """
    n_obs: int: number of observed mutations of a kind
    omega: float: applicable dnds estimate
    omega: float: applicable dnds estimate
    """
    if (n_obs == 0) or np.isnan(n_obs) or np.isnan(omega):
        return n_obs
    elif 0 <= omega <= 1:
        return 0
    elif omega > 1:
        return round((omega - 1) * n_obs / omega)


def excess_rate(n_obs, omega):
    """
    n_obs: int: number of observed mutations of a kind
    omega: float: applicable dnds estimate
    """
    if (n_obs == 0) or np.isnan(n_obs) or np.isnan(omega):
        return 0
    elif 0 <= omega <= 1:
        return 0
    elif omega > 1:
        return (omega - 1) / omega


def add_excess(df):
    df['excess_mis'] = df.apply(lambda v: excess_muts(v['n_mis'], v['wmis_cv']), axis=1)
    df['excess_non'] = df.apply(lambda v: excess_muts(v['n_non'], v['wnon_cv']), axis=1)
    df['excess_spl'] = df.apply(lambda v: excess_muts(v['n_spl'], v['wspl_cv']), axis=1)
    df['excess_rate_mis'] = df.apply(lambda v: excess_rate(v['n_mis'], v['wmis_cv']), axis=1)
    df['excess_rate_non'] = df.apply(lambda v: excess_rate(v['n_non'], v['wnon_cv']), axis=1)
    df['excess_rate_spl'] = df.apply(lambda v: excess_rate(v['n_spl'], v['wspl_cv']), axis=1)
    return df


# def set_role(data, distance_threshold=0.1):
    # """Set the role according to the DNDS output"""
    # if data['wmis_cv'] < 1 and data['wnon_cv'] < 1:  # threshold
        # return "ambiguous"
    # # Check wmis
    # wmis = data['wmis_cv']
    # if wmis >= 1 and data["n_mis"] == 0:
        # wmis = 1
# 
    # # Check wnon
    # wnon = data['wnon_cv']
    # if wnon >= 1 and data["n_non"] == 0:
        # wnon = 1
    # # Those cases with w_non and w_mis <=1 are not informative
    # if wnon <= 1 and wmis <= 1:
        # return "ambiguous"
# 
    # distance = (wmis - wnon) / math.sqrt(2)
    # if distance_threshold is not None and abs(distance) < distance_threshold:
        # return "ambiguous"
    # else:
        # if distance > 0:
            # return 'Act'
        # elif distance < 0:
            # return 'LoF'
        # else:
            # return "ambiguous"

def get_role_cgc (role_cgc):
    if 'oncogene' in role_cgc:
        new_role = 'Act'
    elif 'TSG' in role_cgc:
        new_role = 'LoF'
    else:
        new_role = 'ambiguous'
    return new_role


def set_consensous_role(row):
    if row["ROLE_INTOGEN"] == row["ROLE_CGI"] or row["ROLE_CGI"] == "Unknown":
        return row["ROLE_INTOGEN"]
    else:
        return row["ROLE_CGI"]


def role(df_drivers_role): #, threshold=0.1):

    # # read drivers
    # dndscv_pan_run = bgdata.get('intogen/dndscv/pan')
    # df_drivers_role = pd.read_csv(dndscv_pan_run, sep="\t")
    # df_drivers_role = add_excess(df_drivers_role)
    # df_drivers_role["ROLE_INTOGEN"] = df_drivers_role.apply(lambda row: set_role(row, distance_threshold=threshold),
                                                            # axis=1)

    df_drivers_role = df_drivers_role.rename(columns={"ROLE":"ROLE_INTOGEN"})
    
    # # read mode of action
    # moa = os.path.join(os.environ['INTOGEN_DATASETS'], 'postprocess', 'gene_MoA.tsv')
    # df_moa = pd.read_csv(moa, sep="\t")
    # df_moa.rename(columns={"gene_MoA": "ROLE_CGI"}, inplace=True)

    # read cgc_table, to extract the mode of action / role
    cgc_file = os.path.join(os.environ['INTOGEN_DATASETS'], 'cgc', 'cancer_gene_census_parsed.tsv')
    cgc_df = pd.read_csv(cgc_file,sep='\t')    
    role_df = cgc_df[['Gene Symbol','Role in Cancer']].drop_duplicates()
    role_df['Role in Cancer'].fillna('ambiguous',inplace=True)
    role_df['ROLE_CGI'] = role_df['Role in Cancer'].apply(lambda role_cgc: get_role_cgc(role_cgc) )    

    # Combine both roles
    df_combined_role = pd.merge(df_drivers_role[["SYMBOL", "ROLE_INTOGEN"]], role_df, how="left",
                                left_on=["SYMBOL"], right_on=["Gene Symbol"])
    df_combined_role.drop("Gene Symbol", axis=1, inplace=True)
    df_combined_role.fillna("Unknown", inplace=True)
    df_combined_role["COMBINED_ROLE"] = df_combined_role.apply(lambda row: set_consensous_role(row), axis=1)

    # Update drivers
    df_combined_role.rename(columns={"COMBINED_ROLE": "ROLE"}, inplace=True)
    return df_combined_role[["SYMBOL", "ROLE"]].drop_duplicates()
