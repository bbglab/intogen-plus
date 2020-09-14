
from bgreference import refseq

from intogen_core.formatters.utils import NUCLEOTIDES
from intogen_core.readers import TSVReader

HEADER = ["CHROMOSOME", "POSITION", "REF", "ALT", "SAMPLE", "ID", "GENE", "CONTEXT", "MUTATION_TYPE"]

PYRIMIDINES = {'C', 'T'}
BASE_COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def parse(file):

    for m in TSVReader(file):

        identifier, sample, ref, alt, position = m['#Uploaded_variation'].split('__')
        chromosome, _ = m['Location'].split(':')

        if ref not in NUCLEOTIDES or alt not in NUCLEOTIDES:
            # insertion, deletion or MNV
            continue

        context = refseq('hg38', chromosome, int(position)-1, size=3)
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

        fields = [
            chromosome,
            position,
            ref,
            alt,
            sample,
            identifier,
            m['SYMBOL'],
            context,
            mutation_type
        ]
        yield fields
