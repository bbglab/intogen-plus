import os
import csv
import gzip
import logging

from os import path

LOG = logging.getLogger('intogen')


class VepTask:

    def __init__(self, name, output_folder, conda_env=None, vep_data=None):

        self.name = name
        self.conda_env = conda_env
        self.vep_data = vep_data

        os.makedirs(path.join(output_folder, 'vep'), exist_ok=True)

        self.in_file = path.join(output_folder, "vep", "{}.in.gz".format(name))
        self.in_skip = path.exists(self.in_file)
        self.in_fd = None
        self.in_writer = None

        self.out_file = path.join(output_folder, "vep", "{}.out.gz".format(name))

    def __repr__(self):
        return "Variant Effect Predictor '{}'".format(self.name)

    def input_start(self):
        if not self.in_skip:
            self.in_fd = gzip.open(self.in_file, 'wt')
            self.in_writer = csv.writer(self.in_fd, delimiter='\t')

    def input_write(self, identifier, mut):
        if not self.in_skip:
            # chr - start - end - allele - strand - identifier
            if mut['ALT_TYPE'] == 'snp':
                self.in_writer.writerow( [
                    mut['CHROMOSOME'],
                    mut['POSITION'],
                    mut['POSITION'],
                    "{}/{}".format(mut['REF'], mut['ALT']),
                    mut['STRAND'],
                    identifier
                ])

    def input_end(self):
        if not self.in_skip:
            self.in_writer.close()
            self.in_fd.close()
            self.in_writer = None
            self.in_fd = None

    def run(self):

        # Run vep
        if not path.exists(self.out_file):
            cmd = "bash -c 'source ~/.bashrc && source activate {0} && vep -i <(zcat {1}) -o STDOUT --assembly GRCh37 --no_stats --cache --offline --canonical --dir {3} | gzip > {2}'".format(
                self.conda_env, self.in_file, self.out_file, self.vep_data)
            LOG.debug(cmd)
            os.system(cmd)

        return self.out_file