import os
import sys
import csv
import gzip
import signal
import subprocess

from os import path
from .base import Task


class HotmapsTask(Task):

    KEY = 'hotmapssignature'

    def __init__(self, output_folder, config):

        super().__init__(output_folder, config)

        self.name = None
        self.conda_env = config[self.KEY]['conda_env']
        self.method_folder = config[self.KEY]['method_folder']

        self.in_fd = None
        self.in_writer = None
        self.in_file = None
        self.in_skip = False

        self.out_file = None
        self.output_folder = path.join(output_folder, self.KEY)
        os.makedirs(self.output_folder, exist_ok=True)

    def input_start(self):

        if not self.in_skip:
            self.in_fd = gzip.open(self.in_file, 'wt')
            self.in_writer = csv.writer(self.in_fd, delimiter='\t')
            self.in_writer.writerow(["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele",
                                     "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Variant_Classification", "Transcript_ID", "HGVSp_Short"])

    def input_write(self, _, value):

        if not self.in_skip:

            # Hugo_Symbol Chromosome Start_Position End_Position Reference_Allele Tumor_Seq_Allele2 Tumor_Sample_Barcode
            # Variant_Classification Transcript_ID
            if value['CANONICAL'] == 'YES':
                identifier, sample, ref, alt = value['#Uploaded_variation'].split('__')
                chromosome, position = value['Location'].split(':')
                consequence = value['Consequence'].split(',')[0].replace('missense_variant', 'Missense_Mutation')
                if consequence == "Missense_Mutation":
                    try:
                        aa = value["Amino_acids"].split("/")
                        hgv = "p.{}{}{}".format(aa[0], value['Protein_position'], aa[1])
                    except:
                        hgv = "."
                else:
                    hgv = "."

                self.in_writer.writerow([
                    value['SYMBOL'],
                    chromosome,
                    position,
                    position,
                    ref,
                    alt,
                    sample,
                    consequence,
                    value['Feature'],
                    hgv
                ])

    def input_end(self):
        if not self.in_skip:
            self.in_fd.close()
            self.in_writer = None
            self.in_fd = None

    def run(self):

        # Run HotMaps Signature
        if not path.exists(self.out_file):
            cmd = "bash -c 'source ~/.bashrc && source activate {0} && source {1}/hotmaps.config && {1}/hotmaps.sh {2} {3} {4}'".format(
                self.conda_env, self.method_folder, self.in_file, self.output_folder, os.environ['PROCESS_CPUS'])

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
