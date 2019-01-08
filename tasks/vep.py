import os
import sys
import csv
import gzip
import signal
import subprocess

from os import path
from .base import Task


class VepTask(Task):

    KEY = 'vep'

    def __init__(self, output_folder):

        super().__init__(output_folder)

        self.name = None
        self.vep_data = os.path.join(os.environ['INTOGEN_DATASETS'], 'vep')

        self.in_fd = None
        self.in_writer = None
        self.in_file = None
        self.in_skip = False

        self.out_file = None
        self.output_folder = path.join(output_folder, 'vep')
        os.makedirs(self.output_folder, exist_ok=True)

    def __repr__(self):
        return "Variant Effect Predictor '{}'".format(self.name)

    def input_start(self):
        if not self.in_skip:
            self.in_fd = gzip.open(self.in_file, 'wt')
            self.in_writer = csv.writer(self.in_fd, delimiter='\t')

    def input_write(self, identifier, mut):
        if not self.in_skip:
            # chr - start - end - allele - strand - identifier
            self.in_writer.writerow( [
                mut['CHROMOSOME'],
                mut['POSITION'],
                mut['POSITION'],
                "{}/{}".format(mut['REF'], mut['ALT']),
                mut['STRAND'],
                "{}__{}__{}__{}".format(identifier, mut['SAMPLE'], mut['REF'], mut['ALT'])
            ])

    def input_end(self):
        if not self.in_skip:
            self.in_fd.close()
            self.in_writer = None
            self.in_fd = None

    def run(self):

        # Run vep
        if not path.exists(self.out_file):
            cmd = "bash -c 'singularity run {4} -i <(zcat {0}) -o STDOUT " \
                  "--assembly GRCh37 --no_stats --cache --offline --symbol --protein --tab --canonical --dir {2} " \
                  "| grep -v ^## | gzip > {1}'; " \
                  "mkdir -p vep/logs; mv STDOUT_warnings.txt vep/logs/{3}.log; echo Done;".format(
                self.in_file, self.out_file, self.vep_data, self.name, os.path.join(os.environ['INTOGEN_METHODS'], "vep", "vep.simg"))

            with subprocess.Popen(cmd, shell=True, stdin=sys.stdin, stderr=sys.stderr) as p:
                with open(self.out_file + ".pid", "wt") as fd:
                    fd.write("{}\n".format(p.pid))
                try:
                    errcode = p.wait()
                except:
                    p.kill()
                    p.wait()
                    raise
            os.unlink(self.out_file + ".pid")

            if errcode != 0:
                raise RuntimeError("{} [error] - code {}".format(self, errcode))

        return self.out_file

    def clean(self):

        if path.exists(self.out_file + ".pid"):
            with open(self.out_file + ".pid", "rt") as fd:
                pid = int(fd.readline().strip())
                os.kill(pid, signal.SIGTERM)
            os.unlink(self.out_file + ".pid")

        if path.exists(self.out_file):
            os.unlink(self.out_file)

        if path.exists(self.in_file):
            os.unlink(self.in_file)
