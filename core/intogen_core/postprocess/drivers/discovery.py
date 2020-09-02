import json
import os

import click
import pandas as pd

from intogen_core.postprocess.drivers.bw_list import check_black_white_lists
from intogen_core.postprocess.drivers.data import significative_domains, \
    clusters_2D, clusters_3D, excess
from intogen_core.postprocess.drivers.filters import filter_samples_by_nmuts, \
    filter_by_expression, filter_by_polymorphism, filter_by_olfactory_receptors
from intogen_core.postprocess.drivers.role import role
from intogen_core.postprocess.drivers.signature import analysis_signatures_gene
from intogen_core.postprocess.drivers.vetting import vet


def get_ratio_indels(df_combined):
    df_numbers = df_combined.groupby(["GENE", "TYPE_MUT"], as_index=False).agg({"POSITION": "count"})
    df_agg = df_numbers.pivot_table(columns=["TYPE_MUT"], values=["POSITION"], index=["GENE"],
                                    fill_value=0.0).reset_index()
    df_agg.columns = df_agg.columns.droplevel()
    if df_numbers[df_numbers["TYPE_MUT"] == "INDEL"].shape[0] == 0:
        df_agg["INDEL"] = 0.0
        df_agg.columns = ["GENE", "SNP", "INDEL"]
    else:
        df_agg.columns = ["GENE", "INDEL", "SNP"]

    df_agg["INDEL/SNP"] = df_agg.apply(lambda row: row["INDEL"] / (row["SNP"] + row["INDEL"]), axis=1)

    return df_agg


def include_literature(df):
    # read
    cancermine_file = os.path.join(os.environ['INTOGEN_DATASETS'], 'postprocess', 'cancermine_sentences.tsv')
    cancermine = pd.read_csv(cancermine_file, sep="\t")
    cancermine_g = cancermine.groupby("gene_normalized", as_index=False).agg({"pmid": "count"})
    cancermine_g.rename(columns={"pmid": "n_papers", "gene_normalized": "GENE"}, inplace=True)
    # Match with genes
    df = df.merge(cancermine_g, how="left")
    df["n_papers"].fillna(0, inplace=True)
    return df


def run(combination, mutations, sig_likelihood,
        cohort, ctype,
        smregions, clustl_clusters, hotmaps, dndscv,
            output, muts=3):

    df = pd.read_csv(mutations, sep="\t")

    # 1. Samples with more than 2 mutations in a gene are likely an artifact
    genes_warning_samples = filter_samples_by_nmuts(df, muts)

    # 2. Analysis of signatures
    df_genes, df_combined = analysis_signatures_gene(sig_likelihood, df)

    # 3. Ratio Ratio indels SNP
    df_agg = get_ratio_indels(df_combined)

    # 4. Combine all information
    df = pd.merge(df_agg, df_genes, how="outer")
    df = pd.merge(df, genes_warning_samples, how="left")
    df.fillna(0.0, inplace=True)

    # Filter by expression
    df = filter_by_expression(df, ctype)

    # Filter by Polymorphism
    df = filter_by_polymorphism(df)

     # Add filter by olfactory receptors
    df = filter_by_olfactory_receptors(df)

    # Add filter by known artifacts
    artifacts_file = os.path.join(os.environ['INTOGEN_DATASETS'], 'postprocess', 'artifacts.json')
    with open(artifacts_file) as f:
        artifacts = json.load(f)
    df["Warning_Artifact"] = df.apply(lambda row: row["GENE"] in artifacts["suspects"], axis=1)
    df["Known_Artifact"] = df.apply(lambda row: row["GENE"] in artifacts["known"], axis=1)

    # Add filter by literature
    df = include_literature(df)

    # Checkpoint
    prepared_file = os.path.join(os.path.dirname(output), 'information_vetting_genes.tsv')
    df.to_csv(prepared_file, sep="\t", index=False)


    # Perform vetting
    df = vet(df, combination, ctype)
    df.rename(
        columns={"QVALUE_stouffer_w": "QVALUE_COMBINATION",
                 "Significant_Bidders": "METHODS",
                 "n_papers": "NUM_PAPERS",
                 "cancer_type": "CANCER_TYPE"}, inplace=True)
    df.columns = map(str.upper, df.columns)

    # Checkpoint
    vet_file = os.path.join(os.path.dirname(output), 'all_drivers.tsv')
    columns = ["SYMBOL", "METHODS", "QVALUE_COMBINATION", "TIER",
               "ROLE", "CGC_GENE", "TIER_CGC", "CGC_CANCER_GENE",
               "SIGNATURE9", "WARNING_EXPRESSION", "WARNING_GERMLINE",
               "SAMPLES_3MUTS", "OR_WARNING",
               "KNOWN_ARTIFACT", "NUM_PAPERS", "FILTER"]
    df[columns].sort_values(["SYMBOL"]).to_csv(vet_file, sep="\t", index=False)

    df = check_black_white_lists(df)

    ensembl = os.path.join(os.environ['INTOGEN_DATASETS'], 'regions', 'cds_biomart.tsv')
    df_biomart = pd.read_csv(ensembl, sep="\t", index_col=False, usecols=[1, 10],
                             names=["SYMBOL", "TRANSCRIPT"],
                             header=None)
    df_biomart.drop_duplicates(inplace=True)
    duplicated = df_biomart[df_biomart.duplicated(subset='SYMBOL', keep=False)]['SYMBOL'].unique()

    symbols = df['SYMBOL'].values
    if any(x in symbols for x in duplicated):
        raise Exception('A CGC symbol appears mapped to 2+ transcripts')

    # Checkpoint
    drivers_file = os.path.join(os.path.dirname(output), f'drivers.tsv')
    df.to_csv(drivers_file, sep="\t", index=False)


    # Add drivers data

    df_drivers = df.rename(columns={'MUTS': 'MUTATIONS'})

    df_drivers = df_drivers[['SAMPLES', 'SYMBOL',
                             'METHODS', 'QVALUE_COMBINATION',
                             'CGC_GENE', 'CGC_CANCER_GENE']]
    df_drivers['METHODS'].fillna('combination', inplace=True)

    dfs = [
        significative_domains(smregions),
        clusters_2D(clustl_clusters),
        clusters_3D(hotmaps),
        excess(dndscv),
        role()
    ]

    for df in dfs:  # expected sig. domains, 2D clusters, 3D clusters and excess
        df_drivers = df_drivers.merge(df, how='left')
    # Compute % of samples per cohort
    df_drivers["COHORT"] = cohort

    columns = ["SYMBOL", "COHORT", "METHODS", "SAMPLES",
               "QVALUE_COMBINATION", "CGC_GENE", "CGC_CANCER_GENE",
               "DOMAIN", "2D_CLUSTERS", "3D_CLUSTERS",
               "EXCESS_MIS", "EXCESS_NON", "EXCESS_SPL", "ROLE"]
    df_drivers[columns].sort_values(["SYMBOL"]).to_csv(output, sep="\t", index=False)
    # FIXME TRANSCRIPT


@click.command()
@click.option('--combination', type=click.Path(exists=True), required=True)
@click.option('--mutations', type=click.Path(exists=True), required=True)
@click.option('--sig_likelihood', type=click.Path(exists=True), required=True)
@click.option('--smregions', type=click.Path(exists=True), required=True)
@click.option('--clustl_clusters', type=click.Path(exists=True), required=True)
@click.option('--hotmaps', type=click.Path(exists=True), required=True)
@click.option('--dndscv', type=click.Path(exists=True), required=True)
@click.option('--ctype', type=str, required=True)
@click.option('--cohort', type=str, required=True)
@click.option('--output', type=click.Path(), required=True)
def cli(**kwargs):
    run(**kwargs)


if __name__ == '__main__':
    cli()
