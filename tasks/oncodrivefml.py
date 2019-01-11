import os
import csv
import sys
import gzip
import shutil
import subprocess

from os import path
from .base import Task


class OncodriveFmlTask(Task):

    KEY = 'oncodrivefml'

    def __init__(self, output_folder):
        super().__init__(output_folder)

        self.name = None
        self.in_file = None
        self.in_skip = False
        self.in_fd = None
        self.in_writer = None

        self.out_file = None
        self.output_folder = path.join(output_folder, "oncodrivefml")
        os.makedirs(self.output_folder, exist_ok=True)

    def input_start(self):
        if not self.in_skip:
            self.in_fd = gzip.open(self.in_file, 'wt')
            self.in_writer = csv.writer(self.in_fd, delimiter='\t')
            self.in_writer.writerow(['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE', 'ID'])

    def input_write(self, identifier, mut):
        if not self.in_skip:
            self.in_writer.writerow([
                mut['CHROMOSOME'],
                mut['POSITION'],
                mut['REF'],
                mut['ALT'],
                mut['SAMPLE'],
                identifier
            ])

    def input_end(self):
        if not self.in_skip:
            self.in_fd.close()
            self.in_writer = None
            self.in_fd = None

    def run(self):

        # Run OncodriveFML
        if not path.exists(self.out_file):

            tmp_folder = self.out_file + ".tmp"
            os.makedirs(tmp_folder, exist_ok=True)
           
            cmd = "singularity run {0}/oncodrivefml.simg -i {1} -e {2} -t coding -s {3} -o {4} -c {5} --cores {6}".format(
                os.path.join(os.environ['INTOGEN_METHODS'], 'oncodrivefml'),
                self.in_file,
                os.path.join(os.environ['INTOGEN_DATASETS'], 'oncodrivefml', 'regions_vep88_grch37.gz'),
                'wgs' if "_WGS_" in os.path.basename(self.in_file) else 'wes',
                tmp_folder,
                os.path.join(os.environ['INTOGEN_DATASETS'], 'oncodrivefml', 'oncodrivefml.conf'),
                os.environ.get("INTOGEN_CPUS", 4)
            )

            try:
                o = subprocess.check_output(cmd, shell=True)
            except subprocess.CalledProcessError as e:
                print(e.output.decode())
                sys.exit(e.returncode)
            print(o.decode())

            tmp_output_file = os.path.join(tmp_folder, "{}-oncodrivefml.tsv".format(os.path.basename(self.in_file).split('.')[0]))

            with open(tmp_output_file, "rb") as f_in, gzip.open(self.out_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            shutil.rmtree(tmp_folder)

        return self.out_file

    def clean(self):

        if path.exists(self.out_file):
            os.unlink(self.out_file)

        if path.exists(self.in_file):
            os.unlink(self.in_file)

        tmp_folder = self.out_file + ".tmp"
        if path.exists(tmp_folder):
            shutil.rmtree(tmp_folder)

