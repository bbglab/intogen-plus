import os

import click
import numpy as np
import pandas as pd


def load_cohorts(cohorts):
    df = pd.read_csv(cohorts, sep='\t')
    df.rename(columns={"SAMPLES": "SAMPLES_COHORT", "MUTATIONS": "MUTATIONS_COHORT"}, inplace=True)
    return df[['COHORT', 'CANCER_TYPE', 'SAMPLES_COHORT']]


def load_mutations(mutations):
    df = pd.read_csv(mutations, sep='\t', dtype={'CHR': str})

    mut_counts = df.groupby(['TRANSCRIPT', 'SYMBOL', 'COHORT'], as_index=False).agg({
        'SAMPLES': np.sum
    })
    mut_counts.rename(columns={'SAMPLES': 'MUTATIONS'}, inplace=True)
    return mut_counts


def run(mutations, cohorts, files):
    l = []
    for file in files:
        df = pd.read_csv(file, sep="\t")
        l.append(df)
    drivers = pd.concat(l, sort=True)

    cohort = load_cohorts(cohorts)

    muts = load_mutations(mutations)

    df = pd.merge(drivers, cohort, how='left')
    df = pd.merge(df, muts, how='left')

    # Compute % of samples per cohort
    df["%_SAMPLES_COHORT"] = df.apply(lambda row: row["SAMPLES"] / row["SAMPLES_COHORT"], axis=1)
    # check if the number of mutations is np.nan then discard it
    df = df[~np.isnan(df["MUTATIONS"])]
    df.drop
    df.sort_values(["SYMBOL", "CANCER_TYPE"]).to_csv('drivers.tsv', sep="\t", index=False)

    # Unique drivers
    drivers = df["SYMBOL"].unique()
    # Add the ensembl gene id
    ensembl = os.path.join(os.environ['INTOGEN_DATASETS'], 'regions', 'cds_biomart.tsv')
    df_ensembl = pd.read_csv(ensembl, sep="\t", index_col=False, usecols=[0, 1, 2, 10],
                             names=["ENSEMBL_GENE", "SYMBOL", "ENSEMBL_PROTEIN", "ENSEMBL_TRANSCRIPT"],
                             header=None)  # ENSG00000160752	FDPS	ENSP00000349078	1	155312255	155312395	340	480	1260	1	ENST00000356657
    df_drivers_unique = df_ensembl[df_ensembl["SYMBOL"].isin(drivers)].drop_duplicates()
    df_drivers_unique.to_csv('unique_drivers.tsv', sep="\t", index=False)


@click.command()
@click.option('--mutations', type=click.Path(exists=True), required=True)
@click.option('--cohorts', type=click.Path(exists=True), required=True)
@click.argument('files', nargs=-1)
def cli(mutations, cohorts, files):
    run(mutations, cohorts, files)


if __name__ == "__main__":
    cli()
