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

    def input_write(self, identifier, value):

        if valid_consequence(value['Consequence']) and value['Feature'] in self.transcripts:
            _, sample, ref, alt = value['#Uploaded_variation'].split('__')
            chromosome, position = value['Location'].split(':')

            self.write_row([
                chromosome,
                position,
                ref,
                alt,
                sample,
                self.name,
                identifier
            ])

    def run(self):

        cds_regions = os.path.join(os.environ.get("INTOGEN_DATASETS"), 'shared', 'cds.regions.gz')
        interes_regions = os.path.join(os.environ['INTOGEN_DATASETS'], 'smregions', 'regions_pfam.tsv.gz')
        signatures = '' if self.signatures_file is None else f'-s {self.signatures_file}'

        run_command(
            f"{self.cmdline} -m {self.in_file} -e {cds_regions} -r {interes_regions} {signatures} -o {self.out_file}"
        )

