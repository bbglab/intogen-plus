
import os

from pyliftover import LiftOver

from intogen_core.formatters.utils import NUCLEOTIDES
from intogen_core.readers import TSVReader

HEADER = ["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
          "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode",
          "Variant_Classification", "Transcript_ID", "HGVSp_Short"]

LIFTOVER = LiftOver('hg38', to_db='hg19', search_dir=os.environ['INTOGEN_DATASETS']+'/liftover')


def parse(file):

    for m in TSVReader(file):

        _, sample, ref, alt, position = m['#Uploaded_variation'].split('__')
        chromosome, _ = m['Location'].split(':')

        if ref not in NUCLEOTIDES or alt not in NUCLEOTIDES:
            # insertion, deletion or MNV
            continue

        consequence = m['Consequence'].split(',')[0].replace('missense_variant', 'Missense_Mutation')
        if consequence == "Missense_Mutation":
            try:
                aa = m["Amino_acids"].split("/")
                hgv = "p.{}{}{}".format(aa[0], m['Protein_position'], aa[1])
            except:
                hgv = "."
        else:
            hgv = "."

        strand = '-' if m['STRAND'] == '-1' else '+'
        hg19_position = LIFTOVER.convert_coordinate("chr{}".format(chromosome), int(position) - 1, strand)
        if hg19_position is None or len(hg19_position) != 1:
            continue
        position = hg19_position[0][1] + 1

        fields = [
            m['SYMBOL'],
            chromosome,
            position,
            position,
            ref,
            alt,
            sample,
            consequence,
            m['Feature'],
            hgv
        ]
        yield fields
