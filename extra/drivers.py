import math
import os
import sys
from os import path

import pandas as pd
from bgvep.readers import Tabix


FOLDER = path.dirname(path.abspath(__file__))



def load_chromosomes_genes():
    regions_file = os.path.join(os.environ['INTOGEN_DATASETS'], 'regions', 'cds.regions.gz')
    df = pd.read_csv(regions_file, sep="\t")
    return df[["SYMBOL", "CHROMOSOME", "ELEMENT"]].drop_duplicates()


GENOME_SEQUENCE_MAPS = {'chr{}'.format(c): '{}'.format(c) for c in range(1, 23)}
GENOME_SEQUENCE_MAPS.update({'chrX': 'X', '23': 'X', 'chr23': 'X', 'chrY': 'Y', '24': 'Y', 'chr24': 'Y'})
GENOME_SEQUENCE_MAPS.update({'chrM': 'M', 'MT': 'M', 'chrMT': 'M'})


class TabixAAReader(Tabix):

    def get(self, chromosome, pos, gene):
        chr_ = GENOME_SEQUENCE_MAPS.get(chromosome, chromosome)
        for row in super().get("{}".format(chr_), pos, pos):
            canonical_vep = ((row[-2] == 'YES') and (row[4] == gene))
            if canonical_vep:
                return row[10]


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


def run(output, drivers, vep_version, genome_version):

    df_drivers = pd.read_csv(drivers, sep='\t')
    df_drivers['2D_CLUSTERS'].fillna('', inplace=True)

    chr_df = load_chromosomes_genes()
    df_drivers = df_drivers.merge(chr_df, how='left')

    with TabixAAReader(genome_version, vep_version) as reader:
        df_drivers["2D_CLUSTERS"] = df_drivers.apply(aa_pos, axis=1, reader=reader)

    columns = ["SYMBOL", "TRANSCRIPT", "COHORT", "CANCER_TYPE", "METHODS",
               "MUTATIONS", "SAMPLES", "%_SAMPLES_COHORT",
               "QVALUE_COMBINATION", "ROLE", "CGC_GENE", "CGC_CANCER_GENE",
               "DOMAIN", "2D_CLUSTERS", "3D_CLUSTERS",
               "EXCESS_MIS", "EXCESS_NON", "EXCESS_SPL"]
    df_drivers[columns].sort_values(["SYMBOL", "CANCER_TYPE"]).to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    output = sys.argv[1]
    drivers = sys.argv[2]
    vep = sys.argv[3]
    genome = sys.argv[4]
    run(output, drivers, vep, genome)
