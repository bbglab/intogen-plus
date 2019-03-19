import os
import sys
import csv
import gzip
import signal
import subprocess

from os import path
from .base import Task, run_command

def convert_genome_reference_name(genref):
    if genref == "hg19":
        return "GRCh37"
    return genref.replace("hg", "GRCh")

class VepTask(Task):

    KEY = 'vep'

    def input_write(self, identifier, mut):
        
        # chr - start - end - allele - strand - identifier
        self.in_writer.writerow( [
            mut['CHROMOSOME'],
            mut['POSITION_{}'.format(os.environ['INTOGEN_GENOME'].upper())],
            mut['POSITION_{}'.format(os.environ['INTOGEN_GENOME'].upper())],
            "{}/{}".format(mut['REF'], mut['ALT']),
            mut['STRAND'],
            "{}__{}__{}__{}".format(identifier, mut['SAMPLE'], mut['REF'], mut['ALT'])
        ])

    def run(self):

        cmdline = os.path.join(os.environ['INTOGEN_CONTAINERS'], "{}.simg".format(os.environ['INTOGEN_VEP']))
        assembly = convert_genome_reference_name(os.environ['INTOGEN_GENOME'])
        
        run_command(f"""
            {cmdline} -i {self.in_file} -o STDOUT --assembly {assembly} --no_stats --cache --offline --symbol --protein --tab --canonical --dir {self.datasets} | grep -v ^## | gzip > {self.out_file}; 
            mkdir -p vep/logs; mv STDOUT_warnings.txt vep/logs/{self.name}.log; 
        """)

