import pandas as pd
import numpy as np
import click
import math



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


def set_role(data, distance_threshold=0.1):
    """Set the role according to the DNDS output"""
    if data['wmis_cv'] < 1 and data['wnon_cv'] < 1:  # threshold
        return "ambiguous"
    # Check wmis
    wmis = data['wmis_cv']
    if wmis >= 1 and data["n_mis"] == 0:
        wmis = 1

    # Check wnon
    wnon = data['wnon_cv']
    if wnon >= 1 and data["n_non"] == 0:
        wnon = 1
    # Those cases with w_non and w_mis <=1 are not informative
    if wnon <= 1 and wmis <= 1:
        return "ambiguous"

    distance = (wmis - wnon) / math.sqrt(2)
    if distance_threshold is not None and abs(distance) < distance_threshold:
        return "ambiguous"
    else:
        if distance > 0:
            return 'Act'
        elif distance < 0:
            return 'LoF'
        else:
            return "ambiguous"

def set_consensous_role(row):
    if row["ROLE_INTOGEN"] == row["ROLE_CGI"] or row["ROLE_CGI"] == "Unknown":
        return row["ROLE_INTOGEN"]
    else:
        return row["ROLE_CGI"]







@click.command()
@click.option('--path_cgi_moa',help= 'path to the MoA dataset from cgi', type=click.Path(),required=True)
@click.option('--path_drivers',help= 'path to the drivers from intOGen', type=click.Path(),required=True)
@click.option('--path_run_dndscv',help= 'path to the pan cancer run of dndscv', type=click.Path(),required=True)
@click.option('--path_output',help= 'output path', type=click.Path(),required=True)
@click.option('--threshold',help= 'threshold to be used to make the classification. Default 0.1.',required=False,default=0.1)
def cmdline(path_cgi_moa, path_drivers, path_run_dndscv, path_output, threshold):

    threshold = threshold
    # read drivers
    df_drivers = pd.read_csv(path_drivers, sep="\t")
    drivers = set(df_drivers["SYMBOL"].values)
    df_with_excess_role = pd.read_csv(path_run_dndscv, sep="\t", compression="gzip")
    df_drivers_role = df_with_excess_role[df_with_excess_role["gene_name"].isin(drivers)]
    df_drivers_role = add_excess(df_drivers_role.copy())
    df_drivers_role["ROLE_INTOGEN"] = df_drivers_role.apply(lambda row: set_role(row, distance_threshold=threshold),
                                                            axis=1)
    # read mode of action
    df_moa = pd.read_csv(path_cgi_moa, sep="\t")
    df_moa.rename(columns={"gene_MoA": "ROLE_CGI"}, inplace=True)
    # Combine both roles
    df_combined_role = pd.merge(df_drivers_role[["gene_name", "ROLE_INTOGEN"]], df_moa, how="left",
                                left_on=["gene_name"], right_on=["gene"])
    df_combined_role.drop("gene",axis=1,inplace=True)
    df_combined_role.fillna("Unknown", inplace=True)
    df_combined_role["COMBINED_ROLE"] = df_combined_role.apply(lambda row: set_consensous_role(row), axis=1)
    # Save it
    df_combined_role.to_csv(path_output,sep="\t",index=False)
    # Update drivers
    df_combined_role.rename(columns={"gene_name":"SYMBOL"},inplace=True)
    df_drivers = df_drivers.merge(df_combined_role[["SYMBOL","COMBINED_ROLE"]].drop_duplicates())
    df_drivers.to_csv(path_drivers,sep="\t",index=False)

if __name__ == "__main__":
    cmdline()






