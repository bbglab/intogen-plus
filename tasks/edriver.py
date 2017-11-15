import os
import sys
import csv
import gzip
import signal
import subprocess

from os import path
from .base import Task


class EDriverTask(Task):

    KEY = 'edriver'

    def __init__(self, output_folder, config):

        super().__init__(output_folder, config)

        self.name = None
        self.in_fd = None
        self.in_writer = None
        self.in_file = None
        self.in_skip = False

        self.out_file = None
        self.output_folder = path.join(output_folder, self.KEY)
        os.makedirs(self.output_folder, exist_ok=True)

    def input_start(self):

        if not self.in_skip:
            self.in_fd = gzip.open(self.in_file, 'wt')
            self.in_writer = csv.writer(self.in_fd, delimiter='\t', lineterminator='\n')

    def input_write(self, _, value):

        if not self.in_skip:

            ensp, position = value['ENSP'], value['Protein_position']
            if ensp != "-" and position != "-":
                identifier, sample, ref, alt = value['#Uploaded_variation'].split('__')

                # ENSP Protein_position sample tissue
                self.in_writer.writerow([
                    ensp,
                    position,
                    sample,
                    "unknown"
                ])

    def input_end(self):
        if not self.in_skip:
            self.in_fd.close()
            self.in_writer = None
            self.in_fd = None

    def run(self):

        # Run vep
        if not path.exists(self.out_file):
            cmd = "mkdir -p {2}/tmp_{3} && " \
                  "zcat {1} > {2}/tmp_{3}/input.txt && " \
                  "singularity run {0}/edriver.simg 1 {2}/tmp_{3}/input.txt {0}/iur_and_pfam_ensembl_v85.txt {0}/ensembl_v85_seqs_parsed_header.txt {2}/tmp_{3}/ed &&" \
                  "cat <(head -n1 {2}/tmp_{3}/ed_unknown_with_corrected_p_values.txt) <(grep unknown {2}/tmp_{3}/ed_unknown_with_corrected_p_values.txt) | cut -f2- | gzip > {2}/{3}.out.gz &&" \
                  "rm -rf {2}/tmp_{3}".format(
                os.environ['EDRIVER_DATA'], self.in_file, os.path.abspath(self.output_folder), self.name
            )

            stdout = subprocess.check_output(cmd, shell=True)
            print(stdout.decode())

        return self.out_file
