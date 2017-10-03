import os
import argparse


def parse_arguments():
    info = 'Formats mutation counts on structure to a format consistent with MuPIT table'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-d', '--data-dir',
                        type=str, default='/home/pipeline/mupit_update/tcga/',
                        help='Directory with mutation info (Default: /home/pipeline/mupit_update/tcga/)')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # Process only desired tissue.
    for filename in os.listdir(opts['data_dir']):
        if filename[:9] == 'collected':
            #print filename
            tissue = filename.split('.')[1]
            file_path = os.path.join(opts['data_dir'], filename)

            # make sure it exists
            if os.path.exists(file_path) == False:
                print '  does not exist.'
                continue

            # read in data
            with open(file_path) as handle:
                tcga_lines = handle.readlines()
            num_pdb = len(tcga_lines)

            # write the reformatted output to a file
            out_path = os.path.join(opts['data_dir'], 'mutation_tcga.'+tissue+'.txt')
            with open(out_path, 'w') as wf_mutation:
                # Adding header
                wf_mutation.write('\t'.join(["structure_id", "tissue", "residues", "occurrence"]) + '\n')
                for tcga_line in tcga_lines:
                    [pdb, no_sample, tcga_residues] = tcga_line[:-1].split('\t')
                    tcga_residues_list = tcga_residues.split(', ')
                    for res_ocur in tcga_residues_list:
                        [res, ocur] = res_ocur.split('_')
                        tmp_line = '\t'.join([pdb, tissue, res, ocur]) + '\n'
                        wf_mutation.write(tmp_line)


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
