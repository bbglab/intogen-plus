
from intogen_core.readers import TSVReader

HEADER = ['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE']


def parse(file):

    for m in TSVReader(file):
        fields = [
            m['CHROMOSOME'],
            m['POSITION'],
            m['REF'],
            m['ALT'],
            m['SAMPLE']
        ]
        yield fields
