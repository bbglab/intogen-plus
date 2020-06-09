
from bgreference import refseq

from intogen_core.readers import TSVReader

HEADER = ["CHROMOSOME", "POSITION", "REF", "ALT", "SAMPLE", "ID", "GENE", "CONTEXT", "MUTATION_TYPE"]

PYRIMIDINES = {'C', 'T'}
BASE_COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def parse(file):

    for m in TSVReader(file):

        identifier, sample, ref, alt = m['#Uploaded_variation'].split('__')
        chromosome, position = m['Location'].split(':')

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
