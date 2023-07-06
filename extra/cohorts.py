
import json
import sys
from os import path, environ

import pandas as pd
from bgoncotree.main import BGOncoTree

def read_info(folder):
    info_file = path.join(folder, "info_datasets.tsv")
    df = pd.read_csv(info_file, sep="\t")
    return df


def run(paths, cohorts, dict_label_names, output='cohorts.tsv'):
    data = []
    for path_ in paths:
        data.append(read_info(path_))
    df_info = pd.concat(data, sort=True)
    df_info.columns = map(str.upper, df_info.columns)

    df = pd.read_csv(cohorts, sep='\t')
    df_final = pd.merge(df, df_info, how='left')

    # read labels
    with open(dict_label_names,'r') as f:
        d = json.load(f)

    tree = BGOncoTree(file=path.join(environ['INTOGEN_DATASETS'], 'oncotree', 'tree.tsv'))
    # tree = BGOncoTree(file='/workspace/projects/intogen_plus/fixdatasets-20220321/intogen-plus/datasets/oncotree/tree_old.tsv')
    d = {node.id: node.name for node in tree.descendants('CANCER')}

    df_final["CANCER_TYPE_NAME"] = df_final["CANCER_TYPE"].map(d)
    columns = ["COHORT", "CANCER_TYPE", "CANCER_TYPE_NAME", "SOURCE", "PLATFORM", "REFERENCE", "TYPE", "TREATED", "AGE",  "SAMPLES", "MUTATIONS"]
    df_final[columns].sort_values("CANCER_TYPE").to_csv(output, index=False, sep='\t')


if __name__ == "__main__":
    output = sys.argv[1]
    cohorts = sys.argv[2]
    dict_label_names = sys.argv[3]
    paths = sys.argv[4:]
    run(paths, cohorts, dict_label_names,  output)
