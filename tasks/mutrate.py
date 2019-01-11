import os
import csv
import gzip
import subprocess

from os import path
from .base import Task, valid_consequence


class MutrateTask(Task):

    KEY = 'mutrate'

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

        with open(os.path.join(os.environ['INTOGEN_DATASETS'], 'selected_ensembl_proteins.tsv')) as fd:
            self.proteins = set([r[2] for r in csv.reader(fd, delimiter='\t')])
    
    def run(self):

        # Run vep
        if not path.exists(self.out_file):

            annotmuts = self.in_file.replace(".out.gz", "_annotmuts.out.gz")
            genemuts = self.in_file.replace(".out.gz", "_genemuts.out.gz")
            dndsout = self.in_file
            iterations = 1000
            cores = os.environ.get('INTOGEN_CPUS', 4)

            cmd = "singularity run {0} {1} {2} {3} {4} {5} {6}".format(
                os.path.join(os.environ['INTOGEN_METHODS'], 'mutrate', 'mutrate.simg'),
                annotmuts,
                genemuts,
                dndsout,
                iterations,                
                os.path.join(os.path.abspath(self.output_folder), self.name),
                cores
            )

            stdout = subprocess.check_output(cmd, shell=True)
            print(stdout.decode())

        return self.out_file
