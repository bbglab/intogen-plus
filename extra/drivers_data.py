import glob
import math
import sys
import tempfile
from os import path

import numpy as np
import pandas as pd
import tabix

FOLDER = path.dirname(path.abspath(__file__))


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


def role(dndscv, threshold=0.1, tmp_folder=None):

    # read drivers
    df_drivers_role = pd.read_csv(dndscv, sep="\t")
    df_drivers_role = add_excess(df_drivers_role)
    df_drivers_role["ROLE_INTOGEN"] = df_drivers_role.apply(lambda row: set_role(row, distance_threshold=threshold),
                                                            axis=1)
    # read mode of action
    moa = path.join(FOLDER, 'data', 'gene_MoA.tsv')
    df_moa = pd.read_csv(moa, sep="\t")
    df_moa.rename(columns={"gene_MoA": "ROLE_CGI"}, inplace=True)
    # Combine both roles
    df_combined_role = pd.merge(df_drivers_role[["gene_name", "ROLE_INTOGEN"]], df_moa, how="left",
                                left_on=["gene_name"], right_on=["gene"])
    df_combined_role.drop("gene", axis=1, inplace=True)
    df_combined_role.fillna("Unknown", inplace=True)
    df_combined_role["COMBINED_ROLE"] = df_combined_role.apply(lambda row: set_consensous_role(row), axis=1)

    if tmp_folder:
        comined_role_file = path.join(tmp_folder, 'drivers_role.tsv')
        df_combined_role.to_csv(comined_role_file, sep="\t", index=False)

    # Update drivers
    df_combined_role.rename(columns={"gene_name": "SYMBOL"}, inplace=True)
    return df_combined_role[["SYMBOL", "COMBINED_ROLE"]].drop_duplicates()


def significative_domains(paths):
    QVALUE_THRESHOLD = 0.1
    domains = []
    for path_ in paths:
        for file in glob.glob(path.join(path_, 'smregions', '*.out.gz')):
            cohort = path.basename(file).split(".")[0]
            df = pd.read_csv(file, sep="\t")
            df.rename(columns={"HUGO_SYMBOL": "SYMBOL"}, inplace=True)
            df["COHORT"] = cohort
            # Select significants
            df = df[(df['Q_VALUE'] < QVALUE_THRESHOLD) & ((df['OBSERVED_REGION'] / df['MEAN_SIMULATED']) > 1)]
            if df.shape[0] > 0:
                df[["ENSEMBL_TRANSCRIPT", "PFAM_ID", "START", "END"]] = df['REGION'].str.split(':', expand=True)
                df["DOMAIN"] = df["PFAM_ID"] + ":" + df["START"] + ":" + df["END"]
                df = df.groupby(["SYMBOL", "COHORT"], as_index=False).agg({"DOMAIN": lambda x: ','.join(set(x))})
                domains.append(df)

    df = pd.concat(domains, sort=True)
    return df


GENOME_SEQUENCE_MAPS = {'chr{}'.format(c): '{}'.format(c) for c in range(1, 23)}
GENOME_SEQUENCE_MAPS.update({'chrX': 'X', '23': 'X', 'chr23': 'X', 'chrY': 'Y', '24': 'Y', 'chr24': 'Y'})
GENOME_SEQUENCE_MAPS.update({'chrM': 'M', 'MT': 'M', 'chrMT': 'M'})


class TabixAAReader:

    def __init__(self, file):
        self.file = file
        self.tb = None

    def __enter__(self):
        self.tb = tabix.open(self.file)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return True

    def get(self, chromosome, pos, gene):
        chr_ = GENOME_SEQUENCE_MAPS.get(chromosome, chromosome)
        for row in self.tb.query("{}".format(chr_), pos, pos):
            canonical_vep = ((row[-4] == 'YES') and (row[4] == gene))
            if canonical_vep:
                return row[10]


def get_position_aa(reader, c, chr_, gene_id):
    d = c.split(",")
    start = d[0]
    end = d[-1]
    start_aa = reader.get(chr_, int(start), gene_id)
    end_aa = reader.get(chr_, int(end), gene_id)
    if start_aa > end_aa:
        return pd.Series([end_aa, start_aa])
    else:
        return pd.Series([start_aa, end_aa])


def clusters_2D(paths, vep_tabix):
    PVALUE_THRESHOLD = 0.05
    clusters = []
    with TabixAAReader(vep_tabix) as reader:
        for path_ in paths:
            for file in glob.glob(path.join(path_, 'oncodriveclustl', '*.clusters.gz')):
                cohort = path.basename(file).split(".")[0]
                df = pd.read_csv(file, sep="\t")
                df["COHORT"] = cohort
                # Select significants
                df = df[(df['P'] < PVALUE_THRESHOLD)]
                if df.shape[0] > 0:
                    # Get the amino acid coordinates
                    df[["AA_START", "AA_END"]] = df.apply(
                        lambda row: get_position_aa(reader, row["COORDINATES"], row["CHROMOSOME"], row["ENSID"]), axis=1)
                    df["2D_CLUSTERS"] = df["AA_START"] + ":" + df["AA_END"]
                    df = df.groupby(["SYMBOL", "COHORT"], as_index=False).agg(
                        {"2D_CLUSTERS": lambda x: ','.join(set(map(str, x)))})
                    clusters.append(df)

    df = pd.concat(clusters, sort=True)
    return df


def clusters_3D(paths):
    PVALUE_THRESHOLD = 0.05
    clusters = []
    for path_ in paths:
        for file in glob.glob(path.join(path_, 'hotmaps', '*.clusters.gz')):
            cohort = path.basename(file).split(".")[0]
            df = pd.read_csv(file, sep="\t")
            df.rename(columns={"HUGO Symbol": "SYMBOL", "CRAVAT Res": "3D_CLUSTERS"}, inplace=True)
            df["COHORT"] = cohort
            # Select significants
            df = df[(df['q-value'] < PVALUE_THRESHOLD)]
            if df.shape[0] > 0:
                # Get the amino acid coordinates
                df = df.groupby(["SYMBOL", "COHORT"], as_index=False).agg(
                    {"3D_CLUSTERS": lambda x: ','.join(set(map(str, x)))})
                clusters.append(df)

    df = pd.concat(clusters, sort=True)
    return df


def excess(paths):
    excess = []
    for path_ in paths:
        for file in glob.glob(path.join(path_, 'dndscv', '*.out.gz')):
            cohort = path.basename(file).split(".")[0]
            df = pd.read_csv(file, sep="\t")
            df.rename(columns={"gene_name": "SYMBOL"}, inplace=True)
            df["COHORT"] = cohort
            if df.shape[0] > 0:
                df['EXCESS_MIS'] = df.apply(lambda v: excess_rate(v['n_mis'], v['wmis_cv']), axis=1)
                df['EXCESS_NON'] = df.apply(lambda v: excess_rate(v['n_non'], v['wnon_cv']), axis=1)
                df['EXCESS_SPL'] = df.apply(lambda v: excess_rate(v['n_spl'], v['wspl_cv']), axis=1)
                df = df[["SYMBOL", "COHORT", "EXCESS_MIS", "EXCESS_NON", "EXCESS_SPL"]].drop_duplicates()
                excess.append(df)
    df = pd.concat(excess)
    return df


def mutations(muts):
    mutations_ = pd.read_csv(muts, sep='\t', dtype={'CHR': str})

    mut_counts = mutations_.groupby(['TRANSCRIPT', 'SYMBOL', 'COHORT'], as_index=False).agg({
        'SAMPLES': np.sum
    })
    mut_counts.rename(columns={'SAMPLES': 'MUTATIONS'}, inplace=True)
    return mut_counts


def run(paths, drivers, dndscv, vep, muts, tmp_folder=None):
    tmp_folder = tmp_folder or tempfile.TemporaryDirectory().name

    df_drivers = pd.read_csv(drivers, sep='\t', dtype={'MUTS': int, 'SAMPLES': int})
    df_drivers = df_drivers.rename(columns={'MUTS': 'MUTATIONS', 'Significant_Bidders': 'METHODS'})
    df_drivers = df_drivers[['COHORT', 'MUTATIONS', 'ROLE', 'SAMPLES', 'SYMBOL', 'METHODS', 'TIER']]
    df_drivers['METHODS'].fillna('combination', inplace=True)

    dfs = [
        role(dndscv, tmp_folder=tmp_folder),
        significative_domains(paths),
        clusters_2D(paths, vep),
        clusters_3D(paths),
        excess(paths),
        mutations(muts)
    ]

    for df in dfs:  # expected sig. domains, 2D clusters, 3D clusters and excess
        df_drivers = df_drivers.merge(df, how='left')

    df_drivers.to_csv(drivers, sep="\t", index=False)


if __name__ == "__main__":
    drivers = sys.argv[1]
    dndscv = sys.argv[2]
    vep = sys.argv[3]
    muts = sys.argv[4]
    paths = sys.argv[5:-1]
    tmp_folder = sys.argv[-1]
    run(paths, drivers, dndscv, vep, muts, tmp_folder)

