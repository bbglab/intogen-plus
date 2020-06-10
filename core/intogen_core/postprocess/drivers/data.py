import numpy as np
import pandas as pd


def significative_domains(smregions):
    QVALUE_THRESHOLD = 0.1
    df = pd.read_csv(smregions, sep="\t")
    df.rename(columns={"HUGO_SYMBOL": "SYMBOL"}, inplace=True)
    # Select significants
    df = df[(df['Q_VALUE'] < QVALUE_THRESHOLD) & ((df['OBSERVED_REGION'] / df['MEAN_SIMULATED']) > 1)]

    if df.shape[0] > 0:
        df[["ENSEMBL_TRANSCRIPT", "PFAM_ID", "START", "END"]] = df['REGION'].str.split(':', expand=True)
        df["DOMAIN"] = df["PFAM_ID"] + ":" + df["START"] + ":" + df["END"]
        df = df.groupby(["SYMBOL"], as_index=False).agg({"DOMAIN": lambda x: ','.join(set(x))})
    else:
        df = pd.DataFrame(columns=["SYMBOL", "DOMAIN"])

    return df


def get_position_aa(c):
    d = c.split(",")
    start = int(d[0])
    end = int(d[-1])
    if start > end:
        return pd.Series([str(end), str(start)])
    else:
        return pd.Series([str(start), str(end)])


def clusters_2D(oncodriveclustl_clusters):
    PVALUE_THRESHOLD = 0.05
    df = pd.read_csv(oncodriveclustl_clusters, sep="\t")
    # Select significants
    df = df[(df['P'] < PVALUE_THRESHOLD)]
    if df.shape[0] > 0:
        # Get the amino acid coordinates
        df[["C_START", "C_END"]] = df.apply(
            lambda row: get_position_aa(row["COORDINATES"]), axis=1)
        df["2D_CLUSTERS"] = df["C_START"] + ":" + df["C_END"]
        df = df.groupby(["SYMBOL"], as_index=False).agg(
            {"2D_CLUSTERS": lambda x: ','.join(set(map(str, x)))})
    else:
        df = pd.DataFrame(columns=["SYMBOL", "2D_CLUSTERS"])
    return df


def clusters_3D(hotmaps_file):
    PVALUE_THRESHOLD = 0.05
    df = pd.read_csv(hotmaps_file, sep="\t")
    df.rename(columns={"HUGO Symbol": "SYMBOL", "CRAVAT Res": "3D_CLUSTERS"}, inplace=True)
    # Select significants
    df = df[(df['q-value'] < PVALUE_THRESHOLD)]
    if df.shape[0] > 0:
        # Get the amino acid coordinates
        df = df.groupby(["SYMBOL"], as_index=False).agg(
            {"3D_CLUSTERS": lambda x: ','.join(set(map(str, x)))})
    else:
        df = pd.DataFrame(columns=["SYMBOL", "3D_CLUSTERS"])
    return df


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


def excess(dndscv_file):
    df = pd.read_csv(dndscv_file, sep="\t")
    df.rename(columns={"gene_name": "SYMBOL"}, inplace=True)
    if df.shape[0] > 0:
        df['EXCESS_MIS'] = df.apply(lambda v: excess_rate(v['n_mis'], v['wmis_cv']), axis=1)
        df['EXCESS_NON'] = df.apply(lambda v: excess_rate(v['n_non'], v['wnon_cv']), axis=1)
        df['EXCESS_SPL'] = df.apply(lambda v: excess_rate(v['n_spl'], v['wspl_cv']), axis=1)
        df = df[["SYMBOL", "EXCESS_MIS", "EXCESS_NON", "EXCESS_SPL"]].drop_duplicates()
    else:
        df = pd.DataFrame(columns=["SYMBOL", "EXCESS_MIS", "EXCESS_NON", "EXCESS_SPL"])

    return df

