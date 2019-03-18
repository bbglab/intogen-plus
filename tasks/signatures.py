import os
import csv
import sys
import gzip
import shutil
import subprocess

from os import path
from .base import Task, run_command


class SignatureTask(Task):

    KEY = 'signature'
    INPUT_HEADER = ['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE', 'ID']
    PLATFORM = "WXS"

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
        if self.PLATFORM == 'WGS':
            prefix = 'wg'
        else:
            prefix = 'cds'
        normalize_file = os.path.join(os.environ.get("INTOGEN_DATASETS"), 'shared', f'{prefix}.counts.gz')
        cores = os.environ.get("INTOGEN_CPUS", 4)
        output = os.path.splitext(os.path.splitext(self.out_file)[0])[0]

        run_command(f'''
            {self.cmdline} {self.in_file} {normalize_file} {cores} {output}.json
        ''')
