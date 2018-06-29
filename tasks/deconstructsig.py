import os
import sys
import csv
import gzip
import signal
import subprocess

from os import path
from .base import Task
from bgreference import hg19

PYRIMIDINES = {'C', 'T'}
BASE_COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


class DeconstructSigTask(Task):

    KEY = 'deconstructsig'

    def __init__(self, output_folder, config):

        super().__init__(output_folder, config)

        self.name = None

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
            self.in_writer.writerow(["CHROMOSOME", "POSITION", "REF", "ALT", "SAMPLE", "ID", "GENE", "CONTEXT", "MUTATION_TYPE"])

    def input_write(self, _, value):

        if not self.in_skip:

            # Hugo_Symbol Chromosome Start_Position End_Position Reference_Allele Tumor_Seq_Allele2 Tumor_Sample_Barcode
            # Variant_Classification Transcript_ID
            if value['CANONICAL'] == 'YES':
                identifier, sample, ref, alt = value['#Uploaded_variation'].split('__')
                chromosome, position = value['Location'].split(':')

                context = hg19(chromosome, int(position)-1, size=3)
                if ref in PYRIMIDINES:
                    mutation_type = "{}[{}>{}]{}".format(
                        context[0],
                        ref,
                        alt,
                        context[-1]
                    )
                else:
                    mutation_type = "{}[{}>{}]{}".format(
                        context[0],
                        BASE_COMPLEMENT.get(ref, ref),
                        BASE_COMPLEMENT.get(alt, alt),
                        context[-1]
                    )

                self.in_writer.writerow([
                    chromosome,
                    position,
                    ref,
                    alt,
                    sample,
                    identifier,
                    value['SYMBOL'],
                    context,
                    mutation_type
                ])

    def input_end(self):
        if not self.in_skip:
            self.in_fd.close()
            self.in_writer = None
            self.in_fd = None

    def run(self):

        # Run HotMaps Signature
        if not path.exists(self.out_file):
            cmd = "bash -c 'source ~/.bashrc && source activate intogen_deconstructsig && python {0}/deconstructsig/assign_signature_to_mutation.py -i {1} -o {2}'".format(
                os.environ['INTOGEN_METHODS'], self.in_file, self.out_file)

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
