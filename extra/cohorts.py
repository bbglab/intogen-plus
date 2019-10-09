
import sys
from os import path

import pandas as pd


def read_files(folder):
    samples_file = path.join(folder, 'count_samples.txt')
    samples = pd.read_csv(samples_file, names=["SAMPLES", "COHORT"], sep=" ")

    variants_file = path.join(folder, 'count_variants.txt')
    variants = pd.read_csv(variants_file, names=["COHORT", "MUTATIONS"], sep="\t")

    info_file = path.join(folder, "info_datasets.csv")
    info = pd.read_csv(info_file, sep="\t")

    df = samples.merge(variants)
    df = df.merge(info)
    return df


def run(paths, output='stats_cohorts.tsv'):
    data = []
    for path_ in paths:
        data.append(read_files(path_))
    df_final = pd.concat(data, sort=True)
    df_final.to_csv(output, index=False, sep='\t')


if __name__ == "__main__":
    output = sys.argv[1]
    paths = sys.argv[2:]
    run(paths, output)

