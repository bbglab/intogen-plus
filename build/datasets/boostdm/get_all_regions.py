import gzip
import click
import pandas as pd
import numpy as np
from functools import partial

HEADER = ["GENE_ID", "SYMBOL", "PROTEIN_ID", "CHROMOSOME", "START", "END", "CDS_START", "CDS_END", "CDS_LENGHT", 
          "STRAND", "TRANSCRIPT_ID", "START_EXON", "END_EXON"]


def get_cds_from_biomart(bm):
    """ Extract the cds regions from the canonical biomart dataset
    
    :param bm: str, path to biomart dataset
    """

    biomart_df =  pd.read_csv(bm, sep="\t", names=HEADER)
    biomart_df["STRAND"].replace({-1: "-", 1:"+"}, inplace=True)

    return biomart_df[["CHROMOSOME", "START", "END", "STRAND", "GENE_ID", "SYMBOL", "TRANSCRIPT_ID"]]

def extract_splice_sites(group, num_bases=5):
    """Extract the splice sites of a gene by retaining 5bp in the intronic regions between cds.
    
    :param group: pandas.Dataframe
    :param num_bases: int, number of bases of each splice site
    """
    group.sort_values(by=['START', 'END'], inplace=True)
    group.drop_duplicates(inplace=True)
    group['SP_START_5'] = group['START'] - np.insert(np.repeat(num_bases, len(group) - 1), 0, 0)
    group['SP_END_5'] = group['START'] - np.insert(np.repeat(1, len(group) - 1), 0, 0)
    group['SP_START_3'] = group['END'] + np.append(np.repeat(1, len(group) - 1), 0)
    group['SP_END_3'] = group['END'] + np.append(np.repeat(num_bases, len(group) - 1), 0)


    return group

def main(bm, outf):
    """Reads the biomart dataset and extract the cds regions. 
    Computes the splice sites. 
    Writes output --> cds + splice sites regions.
    
    :param bm: str, path to biomart file
    :param outf: str, path to write output
    """

    cds_df = get_cds_from_biomart(bm)
    
    cds_by_transcript = cds_df.groupby('TRANSCRIPT_ID')

    splice_sites = cds_by_transcript.apply(partial(extract_splice_sites))
    splice_sites.sort_values(by=['CHROMOSOME', 'START', 'END'], inplace=True)

    with gzip.open(outf + "/cds-5spli.regions.gz", 'wt') as output:
        output.write('\t'.join(['CHROMOSOME', 'START', 'END', 'STRAND', 'GENE_ID', 'TRANSCRIPT_ID', 'SYMBOL']) + '\n')
        for count, line in splice_sites.iterrows():
            output.write( 
                f"{line['CHROMOSOME']}\t{line['START']}\t" + \
                f"{line['END']}\t{line['STRAND']}\t{line['GENE_ID']}\t" + \
                f"{line['TRANSCRIPT_ID']}\t{line['SYMBOL']}\n")
            
            ## Adding splice variants when applicable
            if line['SP_START_5'] != line['SP_END_5']:
                output.write(
                    f"{line['CHROMOSOME']}\t{line['SP_START_5']}\t" + \
                    f"{line['SP_END_5']}\t{line['STRAND']}\t{line['GENE_ID']}\t" + \
                    f"{line['TRANSCRIPT_ID']}\t{line['SYMBOL']}\n"
                )
            if line['SP_START_3'] != line['SP_END_3']:
                output.write(
                    f"{line['CHROMOSOME']}\t{line['SP_START_3']}\t" + \
                    f"{line['SP_END_3']}\t{line['STRAND']}\t{line['GENE_ID']}\t" + \
                    f"{line['TRANSCRIPT_ID']}\t{line['SYMBOL']}\n"
                )
    
    return

@click.command()
@click.option('--biomart_path', help= 'biomart file path', type=click.Path(),required=True)
@click.option('--out',help= 'path to the output', type=click.Path(),required=True)
def cli(biomart_path, out):
    main(biomart_path, out)

if __name__ == "__main__":
    cli()