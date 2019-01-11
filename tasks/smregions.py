import os
import csv
import sys
import gzip
import subprocess

from os import path
from .base import Task, valid_consequence


class SmregionsTask(Task):

    KEY = 'smregions'

    def __init__(self, output_folder):

        super().__init__(output_folder)

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
            self.in_writer = csv.writer(self.in_fd, delimiter='\t')
            self.in_writer.writerow(['chr',	'pos', 'ref', 'alt', 'sample', 'Cancer_Type', 'id'])

    def input_write(self, identifier, mut):
        if not self.in_skip:
            self.in_writer.writerow([
                mut['CHROMOSOME'],
                mut['POSITION'],
                mut['REF'],
                mut['ALT'],
                mut['SAMPLE'],
                self.name,
                identifier
            ])

    def input_end(self):
        if not self.in_skip:
            self.in_fd.close()
            self.in_writer = None
            self.in_fd = None

    def run(self):

        # Run SMRegions
        if not path.exists(self.out_file):
            
            cmd = "singularity run {0} -i {1} -e {2} -o {3} -r {4} -c {5}".format(
                os.path.join(os.environ['INTOGEN_METHODS'], 'smregions', 'smregions.simg'),
                self.in_file,
                os.path.join(os.environ['INTOGEN_DATASETS'], 'smregions', '02_cds.regions.gz'),
                self.out_file,
                os.path.join(os.environ['INTOGEN_DATASETS'], 'smregions', 'regions_pfam.tsv.gz'),
                os.path.join(os.environ['INTOGEN_DATASETS'], 'smregions', 'hg19.conf')
            )

            try:
                o = subprocess.check_output(cmd, shell=True)
            except subprocess.CalledProcessError as e:
                print(e.output.decode())
                sys.exit(e.returncode)
            print(o.decode())

        return self.out_file

