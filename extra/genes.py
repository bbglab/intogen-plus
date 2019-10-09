import logging
import sys
import tempfile
from os import path

import pandas as pd

BIOMART_HEADER = ["GENE", "SYMBOL", "PROTEIN", "CHROMOSOME", "START_EXON",
                  "END_EXON", "CDS_START", "CDS_END", "LENGTH", "STRAND",
                  "TRANSCRIPT"]


def run(biomart, cgc, output='genes.tsv', tmp_folder=None):
    tmp_folder = tmp_folder or tempfile.TemporaryDirectory().name

    df_biomart = pd.read_csv(biomart, sep="\t", names=BIOMART_HEADER, index_col=None)
    df_biomart["PROTEIN_LENGTH"] = df_biomart.apply(lambda row: row["LENGTH"] // 3, axis=1)
    df_biomart = df_biomart[['SYMBOL', 'GENE', 'PROTEIN', 'TRANSCRIPT', 'LENGTH', 'PROTEIN_LENGTH']]
    df_biomart.drop_duplicates(inplace=True)

    duplicated = df_biomart[df_biomart.duplicated(subset='SYMBOL', keep=False)]['SYMBOL'].unique()
    with open(path.join(tmp_folder, 'duplicated_symbols.txt'), 'w') as fd:
        fd.write('\n'.join(duplicated))

    df_cgc = pd.read_csv(cgc, sep='\t')
    df_cgc.rename({'Gene Symbol': 'SYMBOL', 'cancer_type': 'CGC_DRIVER_IN'}, axis='columns', inplace=True)
    df_cgc = df_cgc[['SYMBOL', 'CGC_DRIVER_IN']]

    cgc_symbols = df_cgc['SYMBOL'].values
    if any(x in cgc_symbols for x in duplicated):
        logging.error('A CGC symbol appears mapped to 2+ transcripts')
        quit(-1)

    df_genes = df_biomart.merge(df_cgc, how='left', on='SYMBOL')
    df_genes['CGC_DRIVER_IN'].fillna('', inplace=True)

    df_genes.to_csv(output, sep='\t', index=False)


if __name__ == "__main__":
    output = sys.argv[1]
    biomart = sys.argv[2]
    cgc = sys.argv[3]
    tmp_folder = sys.argv[-1]
    run(biomart, cgc, output, tmp_folder)
