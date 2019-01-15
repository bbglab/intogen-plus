import os
import sys


def make_mutations_table_data (mupit_dir):
    file_path = os.path.join(mupit_dir,  'mysql.mutations.tcga.txt')
    wf = open(file_path, 'w')

    for filename in os.listdir(mupit_dir):
        if filename[:14] == 'mutation_tcga.':
            tissue = filename.split('.')[0]
            f = open(os.path.join(mupit_dir, filename))
            for line in f:
                wf.write(line)
            f.close()

    wf.close()


if len(sys.argv) != 2:
    print 'Usage: python make_table_data.py output_dir'
    sys.exit()

outdir = sys.argv[1]

make_mutations_table_data(outdir)
