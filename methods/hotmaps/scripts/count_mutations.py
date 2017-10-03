import os
import sys
import argparse


def parse_arguments():
    info = 'Count mutations for each PDB residue'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-d', '--data-dir',
                        type=str, default='/home/pipeline/mupit_update/tcga/',
                        help='Directory with mutation info (Default: /home/pipeline/mupit_update/tcga/)')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    seqres_by_pdb_alltissues = {}
    no_sample_alltissues = {}

    for filename in os.listdir(opts['data_dir']):
        if filename[:6] == 'mupit.':
            tissue = filename.split('.')[2]

            # get file object
            file_path = os.path.join(opts['data_dir'], filename)
            f = open(file_path)

            # iterate through each line
            seqres_by_pdb = {}
            no_sample = {}
            for line in f:
                [pdbid, seqres, sample_gene_type] = line.rstrip().split('\t')
                pdb = pdbid[:-2]
                chain = pdbid[-1]
                [sample, gene, so] = sample_gene_type.split(';')

                # only use missense mutations
                if so != 'Missense_Mutation':
                    continue

                # Num sample in each tissue
                if no_sample.has_key(pdb) == False:
                    no_sample[pdb] = {}
                no_sample[pdb][sample] = True

                # Num sample in all tissues
                if no_sample_alltissues.has_key(pdb) == False:
                    no_sample_alltissues[pdb] = {}
                no_sample_alltissues[pdb][sample] = True

                if seqres_by_pdb.has_key(pdb) == False:
                    seqres_by_pdb[pdb] = {}
                seqres_chain = seqres + ':' + chain
                # Increases tissue specific occurrence.
                if seqres_by_pdb[pdb].has_key(seqres_chain) == False:
                    seqres_by_pdb[pdb][seqres_chain] = 0
                seqres_by_pdb[pdb][seqres + ':' + chain] += 1
                # Increases all-tissue occurrence.
                if seqres_by_pdb_alltissues.has_key(pdb) == False:
                    seqres_by_pdb_alltissues[pdb] = {}
                if seqres_by_pdb_alltissues[pdb].has_key(seqres_chain) == False:
                    seqres_by_pdb_alltissues[pdb][seqres_chain] = 0
                seqres_by_pdb_alltissues[pdb][seqres_chain] += 1
            f.close()

            # sort PDBs
            pdbs = seqres_by_pdb.keys()
            pdbs.sort()

            # write to file
            out_filename = 'collected.' + filename[6:]
            out_filepath = os.path.join(opts['data_dir'], out_filename)
            with open(out_filepath, 'w') as wf:
                for pdb in pdbs:
                    seqress = seqres_by_pdb[pdb].keys()
                    seqress.sort()
                    line_info = [pdb, str(len(no_sample[pdb].keys())),
                                 ', '.join([seqres + '_' + str(seqres_by_pdb[pdb][seqres])
                                            for seqres in seqress])]
                    wf.write('\t'.join(line_info) + '\n')


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
