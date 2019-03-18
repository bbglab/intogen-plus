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
    PLATFORM = "WXS"

    def input_write(self, identifier, mut):
        self.PLATFORM = mut['PLATFORM']
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
        cds_regions = os.path.join(os.environ.get("INTOGEN_DATASETS"), 'shared', 'cds.regions.gz')
        cores = os.environ.get("INTOGEN_CPUS", 4)
        seq = 'wgs' if self.PLATFORM == "WGS" else 'wes'
        signatures = '' if self.signatures_file is None else f'--signature {self.signatures_file}'
        
        # Run OncodriveFML
        run_command(f"{self.cmdline} -i {self.in_file} -e {cds_regions} -t coding -s {seq} -o {self.out_file} {signatures} --cores {cores}")
