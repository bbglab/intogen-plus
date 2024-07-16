import pandas as pd
import click
import os
import tabix



def read_vep(vep_f):
    return tabix.open(vep_f)

def get(vep_file, chromosome, start, stop):
    """iterates through the vep.tsv.gz"""

    tb = read_vep(vep_file)

    # chr_ = self.map.get(chromosome, chromosome) # what is this??
    for row in tb.query("{}".format(chromosome), start, stop):
        yield row 

def saturation(vep_path, driver_df, regions_df):
    """
    Generator 
    
    :params vep_path: int, path where vep.tsv.gz is stored
    :params driver_df: pandas.DataFrame 
    :regions_df: pandas.DataFrame, canonical.regions.gz

    :yield: driver and df
    """
    d_list = driver_df['SYMBOL'].unique()

    for driver in d_list:
        reg = regions_df[regions_df.SYMBOL == driver].sort_values("GENE_ID")

        # for transcript
        l_cases = list()
        for element in reg.TRANSCRIPT_ID.unique():
            regions_t = reg[reg.TRANSCRIPT_ID == element].sort_values("START")

            for i, r in regions_t.iterrows():
                start = r["START"]
                end = r["END"]
                chr_ = str(r["CHROMOSOME"])
                for data in get(vep_path, chr_, int(start) - 1, int(end)):
                    if ("-" != data[22]) and (data[18] == driver):
                        # then it is the canonical transcript
                        l_cases.append([x for x in data[:25]])   
            
            df = pd.DataFrame(l_cases, columns=['Chromosome', 'Position', 'Reference', 'Alternate', 'Gene', 
                                                'Feature', 'Feature_type', 'Consequence', 'cDNA_position', 
                                                'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 
                                                'Existing_variation', 'Impact','Distance', 'Strand', 'Flags', 
                                                'Symbol', 'Symbol source', 'HGNC_ID', 'Canonical',
                                                'Mane_select', 'Mane_plus_clinical', 'ENSP']
                                )
            
            yield driver, df


@click.command()
@click.option('--drivers', type=click.Path(exists=True), required=True)
def cli(drivers):

    drivers_df = pd.read_csv(drivers, sep='\t', low_memory=False)

    r_path = os.path.join(os.environ['INTOGEN_DATASETS'], 'boostdm', 'saturation','cds-5spli.regions.gz')
    regions_df = pd.read_csv(r_path, sep='\t', low_memory=False)

    vep_path = os.path.join(os.environ['INTOGEN_DATASETS'], 'vep', 'vep.tsv.gz')

    for driver, df in saturation(vep_path, drivers_df, regions_df):
        df.to_csv(f'{driver}.vep.gz', index=False, compression='gzip', sep='\t')
        continue

if __name__ == "__main__":
    cli()   