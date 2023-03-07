import pandas as pd
import click
import os
import tabix


def drivers_list(driver_df):
    """return a list of unique drivers"""
    
    return driver_df['SYMBOL'].unique()

def read_vep(vep_f):
    return tabix.open(vep_f)

def get(vep_file, chromosome, start, stop):
    """iterates through the vep.tsv.gz"""

    tb = read_vep(vep_file)

    # chr_ = self.map.get(chromosome, chromosome) # what is this??
    for row in tb.query("{}".format(chromosome), start, stop):
        yield row 

def saturation(vep_path, d_list, regions_df):
    """
    Generator 
    
    :params vep_path: path where vep.tsv.gz is stored
    :params d_list: list of drivers 
    :regions_df: dataframe of cds.regions.gz

    :yield: a tuple with driver and df
    """
    d = dict()

    for driver in d_list:
        reg = regions_df[regions_df.SYMBOL == driver].sort_values("ELEMENT")

        # for transcript
        l_cases = list()
        for element in reg.ELEMENT.unique():
            regions_t = reg[reg.ELEMENT == element].sort_values("START")

            for i, r in regions_t.iterrows():
                start = r["START"]
                end = r["END"]
                chr_ = str(r["CHROMOSOME"])
                for data in get(vep_path, chr_, int(start), int(end)):
                    if "YES" == data[21]:
                        # then it is the canonical transcript
                        l_cases.append([x for x in data[:23]])   
            df = pd.DataFrame(l_cases,
                              columns=['Chromosome', 'Position', 'Reference', 'Alternate', 'Gene', 'Feature', 'Feature_type', 'Consequence',
                                   'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 
                                   'Impact','Distance', 'Strand', 'Flags', 'Symbol', 'Symbol source', 'HGNC_ID', 'Canonical', 'ENSP']
                                )
            
            yield (driver, df)


def write_out(data, drivers_df):
    """writing output {gene}.{tumortype}.vep.gz"""

    drivers_df.drop_duplicates(subset=['SYMBOL', 'CANCER_TYPE'], inplace= True)
    dic = drivers_df.groupby('SYMBOL')['CANCER_TYPE'].apply(list).to_dict()

    driver, df = data

    for ttype in dic[driver]:
        df.to_csv(f'{driver}.{ttype}.vep.gz', index=False, compression='gzip', sep='\t')

@click.command()
@click.option('--drivers', type=click.Path(exists=True), required=True)
def cli(drivers):

    drivers_df = pd.read_csv(drivers, sep='\t', low_memory=False)
    drivers_l = drivers_list(drivers_df)

    r_path = os.path.join(os.environ['INTOGEN_DATASETS'], 'regions', 'cds.regions.gz')
    regions_df = pd.read_csv(r_path, sep='\t', low_memory=False)

    vep_path = os.path.join(os.environ['INTOGEN_DATASETS'], 'vep', 'vep.tsv.gz')

    for data in saturation(vep_path, drivers_l, regions_df):
        write_out(data, drivers_df)

if __name__ == "__main__":
    cli()   