from intogen_core.readers import TSVReader

HEADER = ['seqnames','start','end', 'strand', 'refnuc', 'varnuc', 'sampleID', 'cancer']


def parse(file):

    for m in TSVReader(file):

        chromosome = 'chr' + m['CHROMOSOME']
        fields = [
            chromosome,
            m['POSITION'],
            m['POSITION'],
            m['STRAND'],
            m['REF'],
            m['ALT'],
            m['SAMPLE'],
            m['CANCER']
        ]
        yield fields
