import os
import sys
import csv
import gzip
import signal
import subprocess

from os import path
from pyliftover import LiftOver
from .base import Task, valid_consequence, run_command


class HotmapsTask(Task):

    KEY = 'hotmaps'
    INPUT_HEADER = ["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", 
                    "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Variant_Classification", "Transcript_ID", "HGVSp_Short"]

    LIFTOVER = None

    def input_start(self):
        super().input_start()

        genome = os.environ['INTOGEN_GENOME'].lower()
        if genome != "hg19":
            self.LIFTOVER = LiftOver(os.path.join(os.environ['INTOGEN_DATASETS'], 'preprocess', '{}ToHg19.over.chain.gz'.format(genome)))

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

            if self.LIFTOVER:
                strand = '-' if value['STRAND'] == '-1' else '+'
                hg19_position = self.LIFTOVER.convert_coordinate("chr{}".format(chromosome), int(position) - 1, strand)
                if hg19_position is None or len(hg19_position) != 1:
                    return      
                position = str(hg19_position[0][1] + 1)

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
        signatures = 'None' if self.signatures_file is None else self.signatures_file

        run_command(
            f"{self.cmdline} {self.in_file} {self.output_folder} {signatures}"
        )
