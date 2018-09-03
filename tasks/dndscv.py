import os
import sys
import csv
import gzip
import subprocess

from os import path
from .base import Task


class DndsCvTask(Task):

    KEY = 'dndscv'

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
            self.in_writer = csv.writer(self.in_fd, delimiter='\t')
            self.in_writer.writerow(["sampleID", "chr", "pos", "ref", "mut"])

    def input_write(self, _, value):

        if not self.in_skip:

            self.in_writer.writerow([
                value['SAMPLE'],
                value['CHROMOSOME'],
                value['POSITION'],
                value['REF'],
                value['ALT']
            ])

    def input_end(self):
        if not self.in_skip:
            self.in_fd.close()
            self.in_writer = None
            self.in_fd = None

    def run(self):

        if not path.exists(self.out_file):

            cmd = "singularity exec {0}/dndscv.img Rscript {0}/run.R {1} {2} {3} {4}".format(
                os.environ['DNDSCV_DATA'],
                self.in_file,
                self.out_file,
                self.out_file.replace(".out.gz", "_annotmuts.out.gz"),
                self.out_file.replace(".out.gz", "_genemuts.out.gz")
            )

            try:
                o = subprocess.check_output(cmd, shell=True)
            except subprocess.CalledProcessError as e:
                print(e.output.decode())
                sys.exit(e.returncode)

            print(o.decode())

        return self.out_file

    def clean(self):

        if path.exists(self.out_file):
            os.unlink(self.out_file)

        if path.exists(self.in_file):
            os.unlink(self.in_file)
