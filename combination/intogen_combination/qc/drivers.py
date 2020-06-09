# Import modules
import os
from collections import defaultdict

import pandas as pd


NEGATIVE_SET = os.path.join(os.environ['INTOGEN_DATASETS'], 'others', 'negative_gene_set.tsv')
POSITIVE_SET = os.path.join(os.environ['INTOGEN_DATASETS'], 'cgc', 'cancer_gene_census_parsed.tsv')


def get_positive_set():
    """Get known cancer genes"""
    drivers = defaultdict(set)
    df_positives = pd.read_csv(POSITIVE_SET,sep="\t")
    for index,row in df_positives.iterrows():
        symbol = row["Gene Symbol"]
        drivers['PANCANCER'].add(symbol)
        if str(row["cancer_type"]) == "nan":
            continue
        for ttype in row["cancer_type"].split(","):
            drivers[ttype].add(symbol)
    return drivers


def get_negative_set():
    """Read a negative set file.
    :return: dictionary
    """
    results = {}
    with open(NEGATIVE_SET, 'r') as fd:
        for line in fd:
            tumor, genes = line.strip().split('\t')
            results[tumor] = set(genes.split(','))
    return results


CGC_GENES_PER_TUMOR = get_positive_set()
NEGATIVE_GENES_PER_TUMOR = get_negative_set()
