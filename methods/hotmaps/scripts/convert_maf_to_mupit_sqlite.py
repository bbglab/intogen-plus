import sys
import os
import sqlite3
import argparse
import pandas as pd
import numpy as np
import time

import maf_utils as mu

def parse_arguments():
    info = "Re-format MAF information to that expected for Reference Transcript info"
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-m', '--maf',
                        type=str, required=True,
                        help='MAF file')
    parser.add_argument('-t', '--tumor-type',
                        type=str, required=True,
                        help='Name of tumor type')
    parser.add_argument('-n', '--no-stratify',
                        action='store_true', default=False,
                        help='Flag indicating whether hypermutators are stratified by tumor type')
    parser.add_argument('-mh', '--mysql-host',
                        type=str,
                        default='karchin-db01.icm.jhu.edu',
                        help='Host name for mysql')
    parser.add_argument('-mdb', '--mysql-db',
                        type=str,
                        default='mupit_modbase',
                        help='Database name for mupit')
    parser.add_argument('--mysql-user',
                        type=str, required=True,
                        help='MySQL user name')
    parser.add_argument('--mysql-passwd',
                        type=str, required=True,
                        help='MySQL password')
    parser.add_argument('-mt', '--mut-threshold',
                        type=int, default=None,
                        help='Use flat mutation count for hypermutator filter (Default: None)')
    parser.add_argument('-i', '--cov-dir',
                        type=str,
                        default='data',
                        help='Directory to write coverage info to')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output file containing info needed for multiple testing correction script')
    args = parser.parse_args()
    return vars(args)


def fix_samp_id(mystring):
    """Remove all of the extra ID info from TCGA barcodes."""
    if mystring.startswith('TCGA'):
        return mystring[:12]
    else:
        return mystring


def read_maf(path, tumor_type):
    """Reads in MAF file with several processing steps to be compatible
    with Xena.

    Steps include:
    * removing duplicate mutations
    * filtering out hypermutate samples
    * keeping only missense mutations
    * adding an ID column for each variant
    * getting codon positions
    * fixing header names

    """
    # figure out whether there is a comment line
    with open(path) as handle:
        first_line = next(handle)
        skip_rows = 1 if first_line.startswith('#') else 0

    # read in data frame
    df = pd.read_csv(path, sep='\t', skiprows=skip_rows, dtype={'Tumor_Sample_Barcode': 'str'})

    # drop duplicate mutations
    #df['Tumor_Sample_Barcode_short'] = df['Tumor_Sample_Barcode'].str[:12]
    df['Tumor_Sample_Barcode_short'] = df['Tumor_Sample_Barcode'].apply(fix_samp_id)
    dup_cols = ['Tumor_Sample_Barcode_short', 'Hugo_Symbol', 'Chromosome',
                'Start_Position', 'End_Position', 'Reference_Allele',
                'Tumor_Seq_Allele2']
    df = df.drop_duplicates(dup_cols)

    #####################
    # filter hypermutated samples based on definition
    # from Kandoth et al.
    #####################
    if opts['no_stratify']:
        strat_col = None
    else:
        strat_col = 'tumor_type'

    # WARNING: No need to filter hypermutators, intogen is filtering them.

    #hypermut_samps, num_mut_list = mu.detect_hypermutators(opts['maf'],
    #                                                       samp_colname='Tumor_Sample_Barcode',
    #                                                       stratify_col=strat_col,
    #                                                       mut_threshold=opts['mut_threshold'])
    #df = df[~df['Tumor_Sample_Barcode'].isin(hypermut_samps)].copy()

    # keep only missense mutations
    df_mis = df[df['Variant_Classification']=='Missense_Mutation'].copy()

    # make variant ID column
    df_mis['ID'] = range(len(df_mis))
    df_mis['ID'] = tumor_type + df_mis['ID'].astype(str).copy()

    # fill in other variants with na's
    is_not_snv = (df_mis['Reference_Allele']=='-') | (df_mis['Tumor_Seq_Allele2']=='-')
    #df_mis.loc[df_mis['Variant_Type']!='SNP', 'HGVSp_Short'] = np.nan
    df_mis.loc[is_not_snv, 'HGVSp_Short'] = np.nan

    # get the mutation info
    df_mis['Reference Codon Position'] = df_mis['HGVSp_Short'].str[3:-1].copy()

    # fix small number of errors in HGVS syntax
    has_letter = df_mis['Reference Codon Position'].str.contains('[A-Za-z]').fillna(True)
    df_mis.loc[has_letter, 'Reference Codon Position'] = '-1'
    is_empty = df_mis['Reference Codon Position']==''
    df_mis.loc[is_empty, 'Reference Codon Position'] = '-1'

    # add mut info columns
    df_mis['Reference Codon Position'] = df_mis['Reference Codon Position'].astype(int).copy()
    df_mis['Reference AA'] = df_mis['HGVSp_Short'].str[2:3].copy()
    df_mis['Alternate AA'] = df_mis['HGVSp_Short'].str[-1].copy()

    # figure out whether 'chr' needs to be added
    df_mis['Chromosome'] = df_mis['Chromosome'].astype(str).copy()
    num_chr = df_mis['Chromosome'].str.startswith('chr').sum()
    if num_chr == 0:
        df_mis['Chromosome'] = 'chr'+df_mis['Chromosome'].astype(str).copy()

    # rename columns to what is expected
    rename_dict = {'Hugo_Symbol': 'HUGO symbol',
                   'Transcript_ID': 'Reference Transcript',
                   'Start_Position': 'Position',
                   'Tumor_Sample_Barcode': 'Sample ID',
                   'Reference_Allele': 'Reference base(s)',
                   'Tumor_Seq_Allele2': 'Alternate base(s)'}
    df_mis = df_mis.rename(columns=rename_dict)

    return df_mis


def main(opts):
    # get mysql connection
    retries = 5
    while retries > 0:
        try:
            db = sqlite3.connect(os.environ['HOTMAPS_DB'])
            db.execute('pragma cache_size = -3000000;')
            db.execute('pragma journal_mode=wal;')
            db.execute('pragma query_only = ON;')
            break
        except Exception:
            time.sleep(5)
            retries -= 1

    if retries == 0:
        raise RuntimeError("Impossible to connect")

    cursor = db.cursor()

    # read in MAF file
    maf_df = read_maf(opts['maf'], opts['tumor_type'])
    out_cols = ['ID', 'Sample ID', 'HUGO symbol', 'Reference Transcript',
                'Reference AA', 'Alternate AA', 'Chromosome', 'Position', 'Reference Codon Position',
                'Reference base(s)', 'Alternate base(s)', 'Strand']
                #'Start_Position', 'End_Position']

    # coverage info file
    cov_info_file = open(opts['cov_dir'] + "/coverage_info.txt", 'w')
    cov_info_file.write('\t'.join(("Tumour Type", "Total Coverage", "Bioassembly Coverage", "Homology Coverage", "Total Mutation Count")) + '\n')

    output_path = opts['output']
    with open(output_path, 'w') as wf:
        # parse header info
        new_header = ["pdb_id", "chain", "residue", 'Reference Genomic Position'] + out_cols
        wf.write('\t'.join(new_header) + '\n')

        # iterate over each line
        count = 0
        mapped_bio_count = 0
        mapped_homology_count = 0
        tot_mapped_count = 0
        for j, row in maf_df.iterrows():
            chrom = row['Chromosome']
            start = row['Position']  # value of Start_Position column
            gene = row['HUGO symbol']

            # query genome2pdb
            myquery = (
                "SELECT gp.PDBId, gp.seqRes, gp.pos1 || ',' || gp.pos2 || ',' || gp.pos3 as `Reference Genomic Position` "
                "FROM ( "
                    "SELECT PDBId, seqRes, pos1, pos2, pos3 "
                    "FROM Genome2PDB "
                    "WHERE chr='{mychr}' AND pos1={mypos} "
                    "UNION "
                    "SELECT PDBId, seqRes, pos1, pos2, pos3 "
                    "FROM Genome2PDB "
                    "WHERE chr='{mychr}' AND pos2={mypos} "
                    "UNION "
                    "SELECT PDBId, seqRes, pos1, pos2, pos3 "
                    "FROM Genome2PDB "
                    "WHERE chr='{mychr}' AND pos3={mypos} "
                ") gp, PDB_Info pi "
                "WHERE gp.PDBId=pi.pdbId AND pi.hugo='{mygene}' AND pi.modbase_filtered=1;"
            ).format(mychr=chrom, mypos=start, mygene=gene)
            cursor.execute(myquery)

            # iterate over mappings
            counted_homology = False
            counted_bio = False
            counted = False
            for result in cursor.fetchall():
                # get mapping info
                (pdbid, seqres, genomic_pos) = result
                chain = pdbid[len(pdbid)-1]
                pdbid = pdbid[:len(pdbid)-2]

                # update mapping statistics
                if not counted:
                    tot_mapped_count += 1
                    counted = True
                if not counted_homology:
                    if pdbid.startswith("NP") or pdbid.startswith("ENSP"):
                        mapped_homology_count += 1
                        counted_homology = True
                if not counted_bio:
                    if not (pdbid.startswith("NP") or pdbid.startswith("ENSP")):
                        mapped_bio_count += 1
                        counted_bio = True

                # write line
                line_str = '\t'.join([pdbid, chain, seqres, genomic_pos] + map(str, row[out_cols].tolist())) + '\n'
                wf.write(line_str)

            count += 1
        tot_coverage = 0 if count == 0 else float(tot_mapped_count)/count*100

        bio_coverage = 0 if tot_mapped_count == 0 else float(mapped_bio_count)/(tot_mapped_count)*100
        homology_coverage = 0 if tot_mapped_count == 0 else float(mapped_homology_count)/(tot_mapped_count)*100
        cov_info_file.write('\t'.join((opts['tumor_type'], str(tot_coverage), str(bio_coverage), str(homology_coverage), str(count))) + '\n')

    cursor.close()
    db.close()


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
