
import logging
from os import path

import click
import pandas as pd
from tqdm import tqdm


def preprocess_df(df):
    """
    Preprocesses the DataFrame by selecting relevant columns, splitting and converting data types.
    
    Expected columns:
    ## Uploaded_variation   Location  Allele    Gene    Feature Feature_type    Consequence cDNA_position   CDS_position    Protein_position    Amino_acids Codons  Existing_variation  IMPACT  DISTANCE    STRAND  FLAGS   SYMBOL  SYMBOL_SOURCE   HGNC_ID CANONICAL   MANE_SELECT MANE_PLUS_CLINICAL  ENSP
    """

    
    df = df[df["MANE_SELECT"] != "-"] #TODO: Check if needed. Since we take the data from the VEPprocessed (already filter for MANE/Transcript). Further investigation needed
    df[['ID_MUT', 'SAMPLES', 'REF', 'ALT', 'POS']] = df['#Uploaded_variation'].str.split('__', expand=True)     # i.e. I0000000000__00493087-9d9d-40ca-86d5-936f1b951c93__C__A__1787655
    df[['CHR', '_']] = df['Location'].str.split(':', expand=True) # i.e. 1:1787655

    df['POS'] = df['POS'].astype(int)
    df.rename(columns={'Feature': 'TRANSCRIPT', 'Consequence': 'CONSEQUENCE'}, inplace=True)
    
    return df


def count_unique(grp):
    s = set(list(grp))
    return len(s)


def run(output, files):
    """
    Reads cohort files, preprocesses them, and aggregates the data.
    """
    list_dfs = []
    for file in tqdm(files):
        cohort = path.basename(file).split(".")[0]
        try:
            df = pd.read_csv(file, sep="\t", low_memory=False)
        except Exception:
            continue

        df = preprocess_df(df)

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
