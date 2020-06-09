
from os import path
import logging

import click
import pandas as pd


def get_info(variation, mut):
    #I0000000000__00493087-9d9d-40ca-86d5-936f1b951c93__C__A	1:1787655
    id_mut, sample_id, ref, alt = variation.split("__")
    chr_, pos = mut.split(":")
    output = pd.Series([str(chr_), int(pos), ref, alt, sample_id])
    return output


def count_unique(grp):
    s = set(list(grp))
    return len(s)


def run(output, files):
    list_dfs = []
    for file in files:
        cohort = path.basename(file).split(".")[0]
        try:
            df = pd.read_csv(file, sep="\t", low_memory=False)
        except:
            continue

        ##Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	IMPACTDISTANCE	STRAND	FLAGS	SYMBOL	SYMBOL_SOURCE	HGNC_ID	CANONICAL	ENSP
        df = df[(df["CANONICAL"] == "YES")][["#Uploaded_variation", "Location", "Feature", 'Consequence', 'Protein_position', "SYMBOL"]]
        df[["CHR", "POS", "REF", "ALT", "SAMPLES"]] = df.apply(lambda row: get_info(row["#Uploaded_variation"], row["Location"]), axis=1)
        df.rename(columns={'Feature': 'TRANSCRIPT', 'Consequence': 'CONSEQUENCE'}, inplace=True)
        df["COHORT"] = cohort
        df["MUTATION"] = df['CHR'] + ':' + df['POS'].astype(str) + ':' + df['REF'] + '>' + df['ALT']
        df.drop(labels=["#Uploaded_variation", "Location"], axis=1)
        df = df.groupby(["MUTATION", "COHORT", "CONSEQUENCE", "CHR", "POS", "REF",
                         "ALT", "TRANSCRIPT", 'Protein_position', "SYMBOL"], as_index=False).agg({"SAMPLES": "count"})
        list_dfs.append(df.drop_duplicates())

    df_final = pd.concat(list_dfs)
    # check duplicated mutations
    x=df_final.groupby(["MUTATION"],as_index=False).agg({"TRANSCRIPT": count_unique})
    if x[x["TRANSCRIPT"] > 0].shape[0]:
        logging.warning('A mutation appears to be mapped to 2+ transcripts')

    df_final.to_csv(output, sep='\t', index=False)


@click.command()
@click.option('-o', '--output', type=click.Path(), required=True)
@click.argument('files', nargs=-1)
def cli(output, files):
    run(output, files)


if __name__ == "__main__":
    cli()
