# Import modules
import os
from collections import defaultdict


NEGATIVE_SET = os.path.join(os.environ['INTOGEN_DATASETS'], 'combination', 'negative_gene_set.tsv')
POSITIVE_SET = os.path.join(os.environ['INTOGEN_DATASETS'], 'combination', 'CGCMay17_cancer_types_TCGA.tsv')


def get_positive_set():
    """Get known cancer genes"""
    drivers = defaultdict(set)
    with open(POSITIVE_SET, 'r') as fd:
        for line in fd:
            if line.startswith('Gene'):
                continue
            symbol, cancer_types = line.strip().split('\t')
            cancer_types = cancer_types.split(',')
            for cancer_type in cancer_types:
                drivers['PANCANCER'].add(symbol)
                drivers[cancer_type].add(symbol)
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
