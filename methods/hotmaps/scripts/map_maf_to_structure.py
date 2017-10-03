"""Maps mutations obtained from maf to protein structure."""
import sys
import os
import MySQLdb
import argparse
import csv
import re


def parse_arguments():
    info = 'Map mutations onto protein structures'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-d', '--data-dir',
                        type=str, required=True,
                        help='Directory for data containing MAF files')
    parser.add_argument('-m', '--match-regex',
                        type=str, default='^TCGA.+.maf$',
                        help='Regex which only matches the initial mutation files')
    parser.add_argument('--host',
                        type=str, required=True,
                        help='MySQL host')
    parser.add_argument('--db',
                        type=str, default='mupit_modbase',
                        help='MySQL MuPIT database name (default: mupit_modbase)')
    parser.add_argument('--mysql-user',
                        type=str, required=True,
                        help='MySQL user name')
    parser.add_argument('--mysql-passwd',
                        type=str, required=True,
                        help='MySQL password')
    parser.add_argument('-o', '--output-dir',
                        type=str, required=True,
                        help='Directory to output results after mapping to protein structure')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # make mysql connection
    db = MySQLdb.connect(host=os.environ['MYSQL_HOST'],
                         port=int(os.environ['MYSQL_PORT']),
                         user=os.environ['MYSQL_USER'],
                         passwd=os.environ['MYSQL_PASSWD'],
                         db=os.environ['MYSQL_DB'])
    cursor = db.cursor()

    # iterate through each tumor type
    for filename in os.listdir(opts['data_dir']):
        #if filename.startswith('tcga') and filename.endswith('.maf'):
	print filename
        if re.match(opts['match_regex'], filename):
            tissue = filename.split('.')[1]
            distinct_variants = {}  # keep track to avoid duplicate variants
            print('Working on:' + tissue)

            # get read file and write file
            file_path = os.path.join(opts['data_dir'], filename)
            write_path = os.path.join(opts['output_dir'], 'non_filtered_mupit.'+filename)
	    print (write_path)
            with open(file_path) as f, open(write_path, 'w') as wf:
                myreader = csv.reader(f, delimiter='\t')

                # parse header
                header = next(myreader)
                if header[0].startswith('#'):
                    # skip comment
                    header = next(myreader)
                gene_ix = header.index('Hugo_Symbol')
                chrom_ix = header.index('Chromosome')
                start_ix = header.index('Start_Position')
                end_ix = header.index('End_Position')
                ref_ix = header.index('Reference_Allele')
                alt_ix = header.index('Tumor_Seq_Allele2')
                samp_ix = header.index('Tumor_Sample_Barcode')
                var_class_ix = header.index('Variant_Classification')

                # iterate over each line
                for line in myreader:
                    # skip if not missense mutations
                    effect = line[var_class_ix]
                    if effect != 'Missense_Mutation':
                        continue

                    # parse the line
                    gene, chrom, start, end = line[gene_ix], line[chrom_ix], line[start_ix], line[end_ix]
                    ref_base, alt_base, sample  = line[ref_ix], line[alt_ix], line[samp_ix]

                    # shorten the sample ID if it is a TCGA barcode
                    if sample.startswith('TCGA'):
                        sample = sample[:12]

                    # check if mutation is a duplicate
                    variant_id = sample + '_' + chrom + '_' + start + '_' + end + '_' + ref_base + '_' + alt_base
                    if distinct_variants.has_key(variant_id):
                        #print(filename + ':' + variant_id + ': duplicate')
                        continue
                    else:
                        distinct_variants[variant_id] = True

                    # query for getting mapped structures from a variant
                    tmp_chr_str = str(chrom) if str(chrom).startswith('chr') else 'chr'+str(chrom)
                    myquery = (
                        "SELECT gp.PDBId, gp.seqRes "
                        "FROM ( "
                            "SELECT PDBId, seqRes "
                            "FROM Genome2PDB "
                            "WHERE chr='{mychr}' AND pos1={mypos} "
                            "UNION "
                            "SELECT PDBId, seqRes "
                            "FROM Genome2PDB "
                            "WHERE chr='{mychr}' AND pos2={mypos} "
                            "UNION "
                            "SELECT PDBId, seqRes "
                            "FROM Genome2PDB "
                            "WHERE chr='{mychr}' AND pos3={mypos} "
                        ") gp, PDB_Info pi "
                        "WHERE gp.PDBId=pi.pdbId AND pi.hugo='{mygene}' AND pi.modbase_filtered=1;"
                    ).format(mychr=tmp_chr_str, mypos=start, mygene=gene)
                    cursor.execute(myquery)

                    # iterate through all mappings
                    for result in cursor.fetchall():
                        (pdbid, seqres) = result
                        wf.write('\t'.join([pdbid, seqres, sample+';'+gene+';'+effect]) + '\n')


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
