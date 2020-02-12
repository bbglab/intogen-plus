
import sys
from os import path
import json
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


def run(paths, dict_label_names, output='cohorts.tsv'):
    data = []
    for path_ in paths:
        data.append(read_files(path_))
    df_final = pd.concat(data, sort=True)
    # read labels
    with open(dict_label_names,'r') as f:
        d = json.load(f)
    df_final.columns = map(str.upper, df_final.columns)
    df_final["CANCER_TYPE_NAME"] = df_final.apply(lambda row: d[row["CANCER_TYPE"]], axis=1)
    columns = ["COHORT", "CANCER_TYPE", "CANCER_TYPE_NAME", "SOURCE", "PLATFORM", "PROJECT", "REFERENCE", "TYPE", "TREATED", "AGE",  "SAMPLES", "MUTATIONS", "WEB_SHORT_COHORT_NAME", "WEB_LONG_COHORT_NAME"]
    df_final[columns].sort_values("CANCER_TYPE").to_csv(output, index=False, sep='\t')


if __name__ == "__main__":
    output = sys.argv[1]
    dict_label_names = sys.argv[2]
    paths = sys.argv[3:]
    run(paths, dict_label_names,  output)

