import os
import sys
import csv
import gzip
import signal
import subprocess

from os import path
from .base import Task


class OncodriveClustlTask(Task):

    KEY = 'oncodriveclustl'

    def __init__(self, output_folder):
        super().__init__(output_folder)

        self.name = None

        self.in_file = None
        self.in_skip = False
        self.in_fd = None
        self.in_writer = None

        self.out_file = None
        self.output_folder = path.join(output_folder, OncodriveClustlTask.KEY)
        os.makedirs(self.output_folder, exist_ok=True)

    def input_start(self):
        if not self.in_skip:
            self.in_fd = gzip.open(self.in_file, 'wt')
            self.in_writer = csv.writer(self.in_fd, delimiter='\t')
            self.in_writer.writerow(['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE', 'ID'])

    def input_write(self, identifier, mut):
        if not self.in_skip:
            self.in_writer.writerow([
                mut['CHROMOSOME'],
                mut['POSITION'],
                mut['REF'],
                mut['ALT'],
                mut['SAMPLE'],
                identifier
            ])

    def input_end(self):
        if not self.in_skip:
            self.in_fd.close()
            self.in_writer = None
            self.in_fd = None

    def run(self):

        # Run vep
        if not path.exists(self.out_file):
            cmd = "oncodriveclustl -i {0} -o {1}/{2} -r {3} -c {4} -sim exon_restricted -simw 35 -sw 45 -cw 45 -cmut 2 -emut 2 -kmer {5} --cds --conseq --oncohort -n 10000 &&" \
                  "(cat {1}/{2}/elements_results.txt | gzip > {1}/{2}.out.gz) &&" \
                  "(cat {1}/{2}/clusters_results.tsv | gzip > {1}/{2}_clusters.out.gz) && " \
                  "(cat {1}/{2}/oncohortdrive_results.out | gzip > {1}/{2}_oncohortdrive.out.gz)".format(
                self.in_file,
                self.output_folder,
                self.name,
                os.path.join(os.environ.get("INTOGEN_DATASETS"), 'oncodrivefml', '02_cds.regions.gz'),
                os.environ.get("PROCESS_CPUS", 4),
                (5 if 'SKCM' in self.name else 3)
            )

            try:
                o = subprocess.check_output(cmd, shell=True)
            except subprocess.CalledProcessError as e:
                print(e.output.decode())
                sys.exit(e.returncode)
            print(o.decode())

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
