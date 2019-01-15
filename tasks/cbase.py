import os
import sys
import csv
import gzip
import subprocess

from os import path
from .base import Task, valid_consequence, run_command
from bgreference import refseq


class CBaseTask(Task):

    KEY = 'cbase'

    CBASE_TO_VEP = [

        ("missense", ["missense_variant", "coding_sequence_variant", "conservative_missense_variant", "rare_amino_acid_variant"]),
        ("nonsense", ['stop_gained', 'stop_lost']),
        ("coding-synon", ["synonymous_variant", "stop_retained_variant"]),
        ("intron", ["transcript_amplification", "intron_variant", "INTRAGENIC", "intragenic_variant"]),
        ("utr-3", ['3_prime_UTR_variant']),
        ("utr-5", ['5_prime_UTR_variant', '5_prime_UTR_premature_start_codon_gain_variant']),
        ("IGR", ['TF_binding_site_variant', 'regulatory_region_variant', 'regulatory_region', 'intergenic_variant', 'intergenic_region']),
    ]

    VEP_TO_CBASE = {}
    for m, veps in CBASE_TO_VEP:
        for v in veps:
            VEP_TO_CBASE[v] = m

    CONTEXT = {
        "AAA":  "0", "AAC":  "1", "AAG":  "2", "AAT":  "3",
        "ACA":  "4", "ACC":  "5", "ACG":  "6", "ACT":  "7",
        "AGA":  "8", "AGC":  "9", "AGG": "10", "AGT": "11",
        "ATA": "12", "ATC": "13", "ATG": "14", "ATT": "15",
        "CAA": "16", "CAC": "17", "CAG": "18", "CAT": "19",
        "CCA": "20", "CCC": "21", "CCG": "22", "CCT": "23",
        "CGA": "24", "CGC": "25", "CGG": "26", "CGT": "27",
        "CTA": "28", "CTC": "29", "CTG": "30", "CTT": "31",
        "GAA": "32", "GAC": "33", "GAG": "34", "GAT": "35",
        "GCA": "36", "GCC": "37", "GCG": "38", "GCT": "39",
        "GGA": "40", "GGC": "41", "GGG": "42", "GGT": "43",
        "GTA": "44", "GTC": "45", "GTG": "46", "GTT": "47",
        "TAA": "48", "TAC": "49", "TAG": "50", "TAT": "51",
        "TCA": "52", "TCC": "53", "TCG": "54", "TCT": "55",
        "TGA": "56", "TGC": "57", "TGG": "58", "TGT": "59",
        "TTA": "60", "TTC": "61", "TTG": "62", "TTT": "63"
    }

    INPUT_HEADER = ["Gene", "mut_eff", "mut_nuc", "context"]

    def input_write(self, _, value):
        if value['Feature'] in self.transcripts:

            consequence = value['Consequence'].split(',')[0]
            mut_eff = self.VEP_TO_CBASE.get(consequence, None)

            # Skip other consequences
            if mut_eff is None:
                return

            _, _, _, mut_nuc = value['#Uploaded_variation'].split('__')

            if mut_nuc not in "ACGT":
                return

            chromosome, position = value['Location'].split(':')
            gene = value['SYMBOL']

            ref = refseq(os.environ['INTOGEN_GENOME'], chromosome, int(position) - 1, size=3).upper()
            context = self.CONTEXT.get(ref, None)
            if context is None:
                return

            self.write_row([
                gene,
                mut_eff,
                mut_nuc,
                context
            ])

    def run(self):
    
        run_command(f"""
            mkdir -p {self.out_file}_tmp && \
            zcat {self.in_file} > {self.out_file}_tmp/input.txt && \
            cd {self.out_file}_tmp && \
            {self.cmdline} {self.out_file}_tmp/input.txt {self.datasets} 0 && \
            tail -n+2 {self.out_file}_tmp/q_values_output.txt | gzip > {self.out_file}
        """)


