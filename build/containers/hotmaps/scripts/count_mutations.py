import os
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

                ## Edit. python3 format to build dictionary for num sample x tissue
                # Num sample in each tissue
                if pdb not in no_sample:
                    no_sample[pdb] = no_sample.get(pdb, {sample: True})
                else:
                    no_sample[pdb][sample] = True

                # Num sample in all tissues
                if pdb not in no_sample_alltissues:
                    no_sample_alltissues[pdb] = no_sample_alltissues.get(pdb, {sample: True})
                else:
                    no_sample_alltissues[pdb][sample] = True

                seqres_chain = seqres + ':' + chain

                ## Edit. Initialize dictionary to count per each sample which residue per chain is mutated. 
                # ex. {"pdb_id":{                           {12xy:{
                #           "sample_id":{                       "YTSC-KL-AXYZ":{
                #                   "res_chain": #                        "16:C_1": 17
                #                       }                                       }
                #               }                                  }
                #      }                                    }
                seqres_by_pdb[pdb] = seqres_by_pdb.get(pdb, {sample:{seqres_chain:0}})
            
                ## Edit. Populate dictionary keeping information per sample
                # Increases tissue specific occurrence.
                if sample not in seqres_by_pdb[pdb]:
                    seqres_by_pdb[pdb][sample] = seqres_by_pdb[pdb].get(sample, {seqres_chain:0})
                    if seqres_chain not in seqres_by_pdb[pdb][sample]:
                        seqres_by_pdb[pdb][sample][seqres_chain] = seqres_by_pdb[pdb][sample].get(seqres_chain, 0) + 1
                    else:
                        seqres_by_pdb[pdb][sample][seqres_chain] += 1
                else:
                    if seqres_chain not in seqres_by_pdb[pdb][sample]:
                        seqres_by_pdb[pdb][sample][seqres_chain] = seqres_by_pdb[pdb][sample].get(seqres_chain, 0) + 1
                    else:
                        seqres_by_pdb[pdb][sample][seqres_chain] += 1

                # Increases all-tissue occurrence.
                seqres_by_pdb_alltissues[pdb] = seqres_by_pdb_alltissues.get(pdb, {sample:{seqres_chain:0}})
                if sample not in seqres_by_pdb_alltissues[pdb]:  
                    seqres_by_pdb_alltissues[pdb][sample] = seqres_by_pdb_alltissues[pdb].get(sample, {seqres_chain:0})
                    if seqres_chain not in seqres_by_pdb[pdb][sample]:
                        seqres_by_pdb_alltissues[pdb][sample][seqres_chain] = seqres_by_pdb_alltissues[pdb][sample].get(seqres_chain, 0) + 1
                    else:
                        seqres_by_pdb_alltissues[pdb][sample][seqres_chain] += 1
                else:
                    if seqres_chain not in seqres_by_pdb_alltissues[pdb][sample]:
                        seqres_by_pdb_alltissues[pdb][sample][seqres_chain] = seqres_by_pdb_alltissues[pdb][sample].get(seqres_chain, 0) + 1
                    else:
                        seqres_by_pdb_alltissues[pdb][sample][seqres_chain] += 1
            f.close()

            # sort PDBs
            pdbs = list(seqres_by_pdb.keys())
            pdbs.sort()

            # write to file
            out_filename = 'collected.' + filename[6:]
            out_filepath = os.path.join(opts['data_dir'], out_filename)
            with open(out_filepath, 'w') as wf:
                for pdb in pdbs:
                    ## Edit. iterate per sample so each structure is repeted as many times as many sample associated to the mutation
                    samples = list(seqres_by_pdb[pdb])
                    samples.sort()
                    for sample in samples:
                        seqress = list(seqres_by_pdb[pdb][sample].keys())
                        seqress.sort()
                        line_info = [pdb, 
                                    str(len(no_sample[pdb].keys())),
                                    ', '.join([seqr + '_' + str(seqres_by_pdb[pdb][sample][seqr])
                                            for seqr in seqress])
                                                , sample
                                                ]
                        wf.write('\t'.join(line_info) + '\n')


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
