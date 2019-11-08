import os
import sys
import csv
import gzip
import subprocess

from os import path
from .base import Task, run_command


class DndsCvTask(Task):

    KEY = 'dndscv'    
    INPUT_HEADER = ["sampleID", "chr", "pos", "ref", "mut"]

    def input_write(self, _, value):
        self.write_row([
            value['SAMPLE'],
            value['CHROMOSOME'],
            value['POSITION_{}'.format(os.environ['INTOGEN_GENOME'].upper())],
            value['REF'],
            value['ALT']
        ])

    def run(self):

        out_annotmuts = self.out_file.replace(".out.gz", ".annotmuts.gz")
        out_genemuts = self.out_file.replace(".out.gz", ".genemuts.gz")
        run_command(f"{self.cmdline} {self.in_file} {self.out_file} {out_annotmuts} {out_genemuts}")
