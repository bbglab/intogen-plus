
from collections import defaultdict
from os import path

import numpy as np
import pandas as pd


def assign_type_mut(value):
    # C[C>T]G SNV; T[C>-]A INDEL; T[C>CA]A INDEL
    muts = value[2:-2]
    muts = muts.replace("-", "")
    nts = muts.split(">")
    if len(nts[0]) != len(nts[1]):
        return "INDEL"
    else:
        return "SNV"


def count_percentage_signature(grp):
    signatures_hypermutators = ["Signature.9", "Signature.10"]
    total = len(grp)
    d = defaultdict(int)
    for sig in grp:
        d[sig] += 1
    d_norm = defaultdict(float)
    for s in d.keys():
        d_norm[s] = (d[s] / total)
    return d_norm[signatures_hypermutators[0]], d_norm[signatures_hypermutators[1]]


def analysis_signatures_gene(signature, df):
    # TODO can this file not exits?
    if path.exists(signature):
        df_sigs = pd.read_csv(signature, sep="\t")
        df_sigs.fillna(0.0, inplace=True)
        if 'Unnamed: 0' in df_sigs.columns.values:
            df_sigs.drop(["Unnamed: 0"], axis=1, inplace=True)
        # get the signature with the highest contribution for a particular sample
        df_sigs["signature_max"] = df_sigs.apply(
            lambda row: row.index[list(row.values).index(np.nanmax(row.values[2:]))], axis=1)
        df_combined = pd.merge(df, df_sigs[["Mutation_type", "signature_max", "Sample"]],
                               left_on=["SAMPLE", "MUTATION_TYPE"], right_on=["Sample", "Mutation_type"],
                               how="left")
        df_combined["signature_max"].fillna("", inplace=True)
        df_combined.drop(columns=["Mutation_type", "Sample"], inplace=True)
    else:
        df_combined = df.copy()
        df_combined["signature_max"] = ""
    df_combined["TYPE_MUT"] = df_combined.apply(lambda row: assign_type_mut(row["MUTATION_TYPE"]), axis=1)

    # Assign % of signatures per gene
    df_genes = df_combined[df_combined["TYPE_MUT"] == "SNV"].groupby("GENE", as_index=False).agg(
        {"signature_max": count_percentage_signature, "POSITION": "count"}).sort_values("POSITION",
                                                                                        ascending=False)
    df_genes["Signature9"] = df_genes.apply(lambda row: row["signature_max"][0], axis=1)
    df_genes["Signature10"] = df_genes.apply(lambda row: row["signature_max"][1], axis=1)
    df_genes.drop(columns=["POSITION", "signature_max"], inplace=True)
    return df_genes, df_combined
