import os
import sys
import csv
import gzip
import signal
import subprocess

from os import path
from .base import Task, valid_consequence, run_command


class HotmapsTask(Task):

    KEY = 'hotmapssignature'
    INPUT_HEADER = ["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", 
                    "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Variant_Classification", "Transcript_ID", "HGVSp_Short"]

    def input_write(self, _, value):

        if valid_consequence(value['Consequence']) and value['Feature'] in self.transcripts:
            _, sample, ref, alt = value['#Uploaded_variation'].split('__')
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

            self.write_row([
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

    def run(self):

        run_command( "singularity run {0} {1} {2} {3}".format(
            os.path.join(os.environ['INTOGEN_CONTAINERS'], 'hotmaps.simg'), self.in_file, self.output_folder, os.environ['INTOGEN_CPUS'])
        )

        
