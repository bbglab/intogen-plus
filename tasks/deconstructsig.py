import os
import sys
import csv
import gzip
import signal
import subprocess

from os import path
from .base import Task, valid_consequence, run_command
from bgreference import refseq

PYRIMIDINES = {'C', 'T'}
BASE_COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


class DeconstructSigTask(Task):

    KEY = 'deconstructsig'
    INPUT_HEADER = ["CHROMOSOME", "POSITION", "REF", "ALT", "SAMPLE", "ID", "GENE", "CONTEXT", "MUTATION_TYPE"]

    def input_write(self, _, value):
        
        # Hugo_Symbol Chromosome Start_Position End_Position Reference_Allele Tumor_Seq_Allele2 Tumor_Sample_Barcode
        # Variant_Classification Transcript_ID
        if valid_consequence(value['Consequence']) and value['Feature'] in self.transcripts:
            identifier, sample, ref, alt = value['#Uploaded_variation'].split('__')
            chromosome, position = value['Location'].split(':')

            context = refseq(os.environ['INTOGEN_GENOME'], chromosome, int(position)-1, size=3)
            if ref in PYRIMIDINES:
                mutation_type = "{}[{}>{}]{}".format(
                    context[0],
                    ref,
                    alt,
                    context[-1]
                )
            else:
                mutation_type = "{}[{}>{}]{}".format(
                    BASE_COMPLEMENT.get(context[-1], context[-1]),
                    BASE_COMPLEMENT.get(ref, ref),
                    BASE_COMPLEMENT.get(alt, alt),
                    BASE_COMPLEMENT.get(context[0], context[0])
                )

            self.write_row([
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

    def run(self):
        run_command(f"{self.cmdline} -i {self.in_file} -o {self.out_file}")

        
