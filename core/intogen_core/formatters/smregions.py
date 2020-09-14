from intogen_core.formatters.utils import NUCLEOTIDES
from intogen_core.readers import TSVReader

HEADER = ['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE']


def parse(file):

    for m in TSVReader(file):
        _, sample, ref, alt, position = m['#Uploaded_variation'].split('__')
        chromosome, _ = m['Location'].split(':')

        if ref not in NUCLEOTIDES or alt not in NUCLEOTIDES:
            # insertion, deletion or MNV
            continue

        fields = [
            chromosome,
            position,
            ref,
            alt,
            sample
        ]
        yield fields
