
from intogen_core.readers import TSVReader

HEADER = ["sampleID", "chr", "pos", "ref", "mut"]


def parse(file):

    for m in TSVReader(file):
        fields = [
            m['SAMPLE'],
            m['CHROMOSOME'],
            m['POSITION'],
            m['REF'],
            m['ALT']
        ]
        yield fields

# FIXME can dnds work on any genome?