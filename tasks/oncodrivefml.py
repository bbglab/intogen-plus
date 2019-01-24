import os
import csv
import sys
import gzip
import shutil
import subprocess

from os import path
from .base import Task, run_command


class OncodriveFmlTask(Task):

    KEY = 'oncodrivefml'    
    INPUT_HEADER = ['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE', 'ID']

    def input_write(self, identifier, mut):
        self.write_row([
            mut['CHROMOSOME'],
            mut['POSITION_{}'.format(os.environ['INTOGEN_GENOME'].upper())],
            mut['REF'],
            mut['ALT'],
            mut['SAMPLE'],
            identifier
        ])

    def run(self):

        # Input parameters
        tmp_folder = self.out_file + ".tmp"
        os.makedirs(tmp_folder, exist_ok=True)
        cds_regions = os.path.join(os.environ.get("INTOGEN_DATASETS"), 'shared', 'cds.regions.gz')
        cores = os.environ.get("INTOGEN_CPUS", 4)
        seq = 'wgs' if "_WGS_" in os.path.basename(self.in_file) else 'wes'
        
        # Run OncodriveFML
        run_command(f"{self.cmdline} -i {self.in_file} -e {cds_regions} -t coding -s {seq} -o {tmp_folder} --cores {cores}")

        # Compress output file
        tmp_output_file = os.path.join(tmp_folder, "{}-oncodrivefml.tsv".format(os.path.basename(self.in_file).split('.')[0]))
        with open(tmp_output_file, "rb") as f_in, gzip.open(self.out_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        shutil.rmtree(tmp_folder)

        
