
from intogen_core.readers import TSVReader

HEADER = ['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE']


def parse(file):

    for m in TSVReader(file):
        _, sample, ref, alt = m['#Uploaded_variation'].split('__')
        chromosome, position = m['Location'].split(':')

        fields = [
            chromosome,
            position,
            ref,
            alt,
            sample
        ]
        yield fields
