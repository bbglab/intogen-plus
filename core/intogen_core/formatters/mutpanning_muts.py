import os

from pyliftover import LiftOver

from intogen_core.readers import TSVReader

HEADER = [
        "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
        "Strand", "Variant_Classification", "Variant_Type", "Reference_Allele",
        "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode",
    ]

LIFTOVER = LiftOver('hg38', to_db='hg19', search_dir=os.environ['INTOGEN_DATASETS']+'/liftover')

def parse(file):

    for m in TSVReader(file):
        _, sample, ref, alt, position = m['#Uploaded_variation'].split('__')
        chromosome, _ = m['Location'].split(':')

        start = int(position)

        variant_type = None
        diff = 0
        if ref == '-':
            variant_type = "INS"
            start = start
            diff = 1
        elif alt == '-':
            variant_type = "DEL"
            diff = len(ref) - 1
        elif len(ref) == len(alt) and len(ref) == 1:
            variant_type = "SNP"
        elif len(ref) == len(alt) and len(ref) == 2:
            variant_type = "DNP"
            diff = 1
        elif len(ref) == len(alt) and len(ref) == 3:
            variant_type = "TNP"
            diff = 2
        else:
            continue

        strand = '-' if m['STRAND'] == '-1' else '+'
        hg19_position = LIFTOVER.convert_coordinate("chr{}".format(chromosome), start - 1, strand)
        if hg19_position is None or len(hg19_position) != 1:
            continue
        start = hg19_position[0][1] + 1

        fields = [
            m['SYMBOL'],
            chromosome,
            f"{start}",
            f"{start+diff}",
            m['STRAND'],
            m['Consequence'],
            variant_type,
            ref,
            ref,
            alt,
            sample
        ]
        yield fields
