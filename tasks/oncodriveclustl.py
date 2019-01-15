import os
import sys
import csv
import gzip
import signal
import subprocess

from os import path
from .base import Task, run_command


class OncodriveClustlTask(Task):

    KEY = 'oncodriveclustl'    
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

        cds_regions = os.path.join(os.environ.get("INTOGEN_DATASETS"), 'shared', 'cds.regions.gz')
        cores = os.environ.get("INTOGEN_CPUS", 4)
        genome = os.environ['INTOGEN_GENOME']
        kmer = (5 if 'SKCM' in self.name else 3)

        run_command(f"""
            {self.cmdline} -i {self.in_file} -o {self.output_folder}/{self.name} -r {cds_regions}  
            -c {cores} -g {genome} -sim region_restricted -simw 35 -sw 45 -cw 45 -cmut 2 -emut 2 -kmer {kmer} -n 10000 && 
            (cat {self.output_folder}/{self.name}/elements_results.txt | gzip > {self.out_file}) && 
            (cat {self.output_folder}/{self.name}/clusters_results.tsv | gzip > {self.output_folder}/{self.name}.clusters.gz) &&  
            (cat {self.output_folder}/{self.name}/oncohortdrive_results.out | gzip > {self.output_folder}/{self.name}.oncohortdrive.gz) 
        """)
