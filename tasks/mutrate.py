import os
import csv
import gzip
import subprocess

from os import path
from .base import Task, valid_consequence, run_command


class MutrateTask(Task):

    KEY = 'mutrate'
    
    def run(self):
        
        annotmuts = self.in_file.replace(".out.gz", ".annotmuts.gz")
        genemuts = self.in_file.replace(".out.gz", ".genemuts.gz")
        dndsout = self.in_file
        iterations = 1000
        cores = os.environ.get('INTOGEN_CPUS', 4)

        run_command(f"{self.cmdline} {annotmuts} {genemuts} {dndsout} {iterations} {self.out_file} {cores}")