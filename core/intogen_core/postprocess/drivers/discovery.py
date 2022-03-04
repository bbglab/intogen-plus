import json
import os

import click
import pandas as pd

from intogen_core.exceptions import IntogenError
# from intogen_core.postprocess.drivers.bw_list import check_black_white_lists
#from intogen_core.postprocess.drivers.bw_list import read_file
from intogen_core.postprocess.drivers.data import significative_domains, \
    clusters_2D, clusters_3D, excess
from intogen_core.postprocess.drivers.filters import filter_samples_by_nmuts, \
    filter_by_expression, filter_by_polymorphism, filter_by_olfactory_receptors
from intogen_core.postprocess.drivers.role import role
from intogen_core.postprocess.drivers.signature import analysis_signatures_gene
from intogen_core.postprocess.drivers.vetting import vet


VET_COLUMNS = ["SYMBOL", "ALL_METHODS", "SIG_METHODS", "QVALUE_COMBINATION", "QVALUE_CGC_COMBINATION", 
           "RANKING", "TIER", "ROLE", "CGC_GENE", "TIER_CGC", "CGC_CANCER_GENE",
           "SIGNATURE9", "SIGNATURE10", "WARNING_EXPRESSION", "WARNING_GERMLINE",
           "SAMPLES_3MUTS", "OR_WARNING", "WARNING_ARTIFACT", "KNOWN_ARTIFACT",
           "NUM_PAPERS","WARNING_ENSEMBL_TRANSCRIPTS","DRIVER", "FILTER"]

DRIVERS_COLUMNS = ["SYMBOL", "COHORT", "METHODS", "SAMPLES",
               "QVALUE_COMBINATION", "CGC_GENE", "CGC_CANCER_GENE",
               "DOMAIN", "2D_CLUSTERS", "3D_CLUSTERS",
               "EXCESS_MIS", "EXCESS_NON", "EXCESS_SPL", "ROLE"]


def get_ratio_indels(df_combined):
    df_numbers = df_combined.groupby(["GENE", "TYPE_MUT"], as_index=False).agg({"POSITION": "count"})
    df_agg = df_numbers.pivot_table(columns=["TYPE_MUT"], values=["POSITION"], index=["GENE"],
                                    fill_value=0.0).reset_index()
    df_agg.columns = df_agg.columns.droplevel()
    if df_numbers[df_numbers["TYPE_MUT"] == "INDEL"].shape[0] == 0:
        df_agg["INDEL"] = 0.0
        df_agg.columns = ["GENE", "SNV", "INDEL"]
    else:
        df_agg.columns = ["GENE", "INDEL", "SNV"]

    df_agg["INDEL/SNV"] = df_agg.apply(lambda row: row["INDEL"] / (row["SNV"] + row["INDEL"]), axis=1)

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

def read_file(filein):
    f = open(filein, 'r')
    genes = set()
    for line in f.readlines():
        line = line.strip()
        genes.add(line)
    f.close()
    return genes

def run(combination, mutations, sig_likelihood,
        cohort, ctype,
        smregions, clustl_clusters, hotmaps, dndscv,
            output_drivers, output_vet, muts=3):

    df = pd.read_csv(mutations, sep="\t")

    # 1. Samples with more than 2 mutations in a gene are likely an artifact
    print('1. Samples with more than 2 mutations in a gene are likely an artifact')
    genes_warning_samples = filter_samples_by_nmuts(df, muts)

    # 2. Analysis of signatures: to detect samples with hypermutated behaviour
    print('2. Analysis of signatures: to detect samples with hypermutated behaviour')
    df_genes, df_combined = analysis_signatures_gene(sig_likelihood, df)

    # 3. Ratio Ratio indels SNV
    print('3. Ratio indels SNV')
    df_agg = get_ratio_indels(df_combined)

    # 4. Combine all information
    print('4. Combine all information')
    df = pd.merge(df_agg, df_genes, how="outer")
    df = pd.merge(df, genes_warning_samples, how="left")
    df.fillna(0.0, inplace=True)

    # 5. Add expression info
    print('5. Add expression info')
    df = filter_by_expression(df, ctype)

    # 6. Add Polymorphism info
    print('6. Add Polymorphism info')
    df = filter_by_polymorphism(df)

    # 7. Add filter by olfactory receptor
    print('7. Add filter by olfactory receptors')
    df = filter_by_olfactory_receptors(df)

    # 8. Add filter by known artifacts and black listed
    print('8. Add filter by known artifacts and black listed')
    artifacts_file = os.path.join(os.environ['INTOGEN_DATASETS'], 'postprocess', 'artifacts.json')
    with open(artifacts_file) as f:
        artifacts = json.load(f)
    black_listed_file = os.path.join(os.environ['INTOGEN_DATASETS'], 'postprocess', 'black_listed.txt')
    black_listed = read_file(black_listed_file)

    df["Warning_Artifact"] = df.apply(lambda row: (row["GENE"] in artifacts["suspects"]), axis=1)
    df["Known_Artifact"] = df.apply(lambda row: row["GENE"] in artifacts["known"] or (row['GENE'] in black_listed), axis=1)

    # 9. Add filter by literature (cancermine) and white listed
    print('9. Add filter by literature (cancermine) and white listed')
    df = include_literature(df)

    # 10. Perform vetting
    print('10. Perform vetting')
    # try:
    df = vet(df, combination, ctype)
    # except IntogenError as e:
        # # NO drivers to perform the vetting
        # print(str(e))
        # # Create and empty dataframe and exit without error
        # df = pd.DataFrame(columns=OUT_COLUMNS)
        # df.to_csv(output, sep="\t", index=False)
        # return
    if len(df) == 0:
        vetting_df = pd.DataFrame(columns = VET_COLUMNS)
        vetting_df.to_csv(output_vet, sep="\t", index=False)
        
        drivers_df = pd.DataFrame(columns = DRIVERS_COLUMNS)
        drivers_df.to_csv(output_drivers, sep="\t", index=False)       

    else:
        # 11. Check number of ENSEMBL transcripts per gene
        print('11. Check number of ENSEMBL transcripts per gene')
        ensembl = os.path.join(os.environ['INTOGEN_DATASETS'], 'regions', 'cds_biomart.tsv')
        df_biomart = pd.read_csv(ensembl, sep="\t", index_col=False, usecols=[1, 10],
                                 names=["SYMBOL", "TRANSCRIPT"],
                                 header=None)
        df_biomart.drop_duplicates(inplace=True)
        duplicated = df_biomart[df_biomart.duplicated(subset='SYMBOL', keep=False)]['SYMBOL'].unique()

        symbols = df['SYMBOL'].values
        if any(x in symbols for x in duplicated):
            raise Exception('A CGC symbol appears mapped to 2+ transcripts')

        df['WARNING_ENSEMBL_TRANSCRIPTS'] = df['SYMBOL'].apply(lambda x: x in duplicated)

        # 12. Prepare file with arranged columns
        print('12. Prepare file with arranged columns')
        df.rename(
            columns={"QVALUE_stouffer_w": "QVALUE_COMBINATION",
                     "All_Bidders": "ALL_METHODS",
                     "Significant_Bidders":"SIG_METHODS",
                     "n_papers": "NUM_PAPERS",
                     "cancer_type": "CANCER_TYPE",
                     "MUTS":"MUTATIONS",
                     "QVALUE_CGC_stouffer_w":"QVALUE_CGC_COMBINATION"}, inplace=True)
        df.columns = map(str.upper, df.columns)

        # 13. Checkpoint: save file with vetting info
        print('13. Checkpoint: save file with vetting info')
        
        df['SIG_METHODS'].fillna('combination', inplace=True)
        df[VET_COLUMNS].sort_values(["SYMBOL"]).to_csv(output_vet, sep="\t", index=False)

        # 14. Filter drivers

        df_drivers = df[df['FILTER']=='PASS']
        
        df_drivers = df_drivers.rename(columns={'SIG_METHODS':'METHODS'})

        df_drivers = df_drivers[['SAMPLES', 'SYMBOL',
                                 'METHODS', 'QVALUE_COMBINATION',
                                 'CGC_GENE', 'CGC_CANCER_GENE','ROLE']]

        # 15. Add mutational features
        print('15. Add mutational features')
        dfs = [
            significative_domains(smregions),
            clusters_2D(clustl_clusters),
            clusters_3D(hotmaps),
            excess(dndscv),
            role(df_drivers)
        ]

        for df in dfs:  # expected sig. domains, 2D clusters, 3D clusters and excess
            df_drivers = df_drivers.merge(df, how='left')
        # Compute % of samples per cohort
        df_drivers["COHORT"] = cohort

        df_drivers[DRIVERS_COLUMNS].sort_values(["SYMBOL"]).to_csv(output_drivers, sep="\t", index=False)
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
@click.option('--output_drivers', type=click.Path(), required=True)
@click.option('--output_vet', type=click.Path(), required=True)
def cli(**kwargs):
    run(**kwargs)


if __name__ == '__main__':
    cli()
