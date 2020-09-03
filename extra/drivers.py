import math
import os
import sys
from os import path

import numpy as np
import pandas as pd
import tabix

# from bgvep.readers import Tabix

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


def role(dndscv, threshold=0.1):
    # read drivers
    df_drivers_role = pd.read_csv(dndscv, sep="\t" ,low_memory=False)
    df_drivers_role = add_excess(df_drivers_role)
    df_drivers_role["ROLE_INTOGEN"] = df_drivers_role.apply(lambda row: set_role(row, distance_threshold=threshold),
                                                            axis=1)
    # read mode of action
    moa = path.join(FOLDER, 'data', 'gene_MoA.tsv')
    df_moa = pd.read_csv(moa, sep="\t" ,low_memory=False)
    df_moa.rename(columns={"gene_MoA": "ROLE_CGI"}, inplace=True)
    # Combine both roles
    df_combined_role = pd.merge(df_drivers_role[["gene_name", "ROLE_INTOGEN"]], df_moa, how="left",
                                left_on=["gene_name"], right_on=["gene"])
    df_combined_role.drop("gene", axis=1, inplace=True)
    df_combined_role.fillna("Unknown", inplace=True)
    df_combined_role["COMBINED_ROLE"] = df_combined_role.apply(lambda row: set_consensous_role(row), axis=1)

    # Update drivers
    df_combined_role.rename(columns={"gene_name": "SYMBOL"}, inplace=True)
    return df_combined_role[["SYMBOL", "COMBINED_ROLE"]].drop_duplicates()


def load_chromosomes_genes():
    regions_file = os.path.join(os.environ['INTOGEN_DATASETS'], 'regions', 'cds.regions.gz')
    df = pd.read_csv(regions_file, sep="\t", low_memory=False)
    return df[["SYMBOL", "CHROMOSOME", "ELEMENT"]].drop_duplicates()


GENOME_SEQUENCE_MAPS = {'chr{}'.format(c): '{}'.format(c) for c in range(1, 23)}
GENOME_SEQUENCE_MAPS.update({'chrX': 'X', '23': 'X', 'chr23': 'X', 'chrY': 'Y', '24': 'Y', 'chr24': 'Y'})
GENOME_SEQUENCE_MAPS.update({'chrM': 'M', 'MT': 'M', 'chrMT': 'M'})

class TabixAAReader:

    def __init__(self):
        self.file = vep   # I have put the vep file , it was this: bgdata.get_path('vep', 'wgs_tabix', '{}_{}'.format(genome, vep_build))
        self.tb = None


    def get(self, chromosome, pos, gene):
        chr_ = GENOME_SEQUENCE_MAPS.get(chromosome, chromosome)
        tb = tabix.open(self.file)
        for row in self.tb.query("{}".format(chr_), pos, pos):
            # it was this: for row in super().get("{}".format(chr_), pos, pos):
            if row[4] == gene:
                return row[10]

    def __enter__(self):
        self.tb = tabix.open(self.file)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return False


def get_position_aa(reader, gene_id, chr_, start, end):
    start_aa = reader.get(chr_, start, gene_id)
    end_aa = reader.get(chr_, end, gene_id)
    if int(start_aa) > int(end_aa):
        return ':'.join([end_aa, start_aa])
    else:
        return ':'.join([start_aa, end_aa])


def aa_pos(x, reader):
    if x["2D_CLUSTERS"] == '':
        return ''
    l = []
    for cluster in x["2D_CLUSTERS"].split(','):
        i, j = cluster.split(':')
        start, end = int(i), int(j)
        l.append(get_position_aa(reader, x["ELEMENT"], x["CHROMOSOME"], start, end))

    return ','.join(l)


# def run(output, drivers, dndscv, vep_version, genome_version):
def run(output, drivers):
    df_drivers = pd.read_csv(drivers, sep='\t', low_memory=False)
    df_drivers['2D_CLUSTERS'].fillna('', inplace=True)

    # df_role = role(dndscv)
    # df_drivers = df_drivers.merge(df_role, how='left')

    chr_df = load_chromosomes_genes()
    df_drivers = df_drivers.merge(chr_df, how='left')
    # with TabixAAReader(genome_version, vep_version) as reader:
    with TabixAAReader() as reader:
        df_drivers["2D_CLUSTERS"] = df_drivers.apply(aa_pos, axis=1, reader=reader)

    #df_drivers.rename(columns={"COMBINED_ROLE": "ROLE"}, inplace=True)

    #columns = ["SYMBOL", "TRANSCRIPT", "COHORT", "CANCER_TYPE", "METHODS",
     #          "MUTATIONS", "SAMPLES", "%_SAMPLES_COHORT",
      #         "QVALUE_COMBINATION", "ROLE", "CGC_GENE", "CGC_CANCER_GENE",
       #        "DOMAIN", "2D_CLUSTERS", "3D_CLUSTERS",
        #       "EXCESS_MIS", "EXCESS_NON", "EXCESS_SPL"]
    columns = ["SYMBOL", "TRANSCRIPT", "COHORT", "CANCER_TYPE", "METHODS",
               "MUTATIONS", "SAMPLES", "%_SAMPLES_COHORT",
               "QVALUE_COMBINATION", "CGC_GENE", "CGC_CANCER_GENE",
               "DOMAIN", "2D_CLUSTERS", "3D_CLUSTERS",
               "EXCESS_MIS", "EXCESS_NON", "EXCESS_SPL"]

    df_drivers[columns].sort_values(["SYMBOL", "CANCER_TYPE"]).to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    output = sys.argv[1]
    drivers = sys.argv[2]
    # dndscv = sys.argv[3]
    vep = sys.argv[3]  # not necessary
    #genome = sys.argv[4]  # not necessary
    # run(output, drivers, dndscv, vep, genome)
    run(output, drivers)

# how to run?
# python drivers.py drivers.out drivers.tsv VEP_canonical_transcripts.out.gz 
