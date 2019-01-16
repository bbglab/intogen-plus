import os
import csv
import sys
import gzip
import subprocess

from os import path
from .base import Task, valid_consequence, run_command


class SmregionsTask(Task):

    KEY = 'smregions'
    INPUT_HEADER = ['chr',	'pos', 'ref', 'alt', 'sample', 'Cancer_Type', 'id']

    def input_write(self, identifier, mut):
        self.write_row([
            mut['CHROMOSOME'],
            mut['POSITION_{}'.format(os.environ['INTOGEN_GENOME'].upper())],
            mut['REF'],
            mut['ALT'],
            mut['SAMPLE'],
            self.name,
            identifier
        ])

    def run(self):

        cds_regions = os.path.join(os.environ.get("INTOGEN_DATASETS"), 'shared', 'cds.regions.gz')
        interes_regions = os.path.join(os.environ['INTOGEN_DATASETS'], 'smregions', 'regions_pfam.tsv.gz')
    
        run_command(f"{self.cmdline} -m {self.in_file} -e {cds_regions} -r {interes_regions} -o {self.out_file}")

