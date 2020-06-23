import os

from intogen_core.readers import TSVReader

HEADER = [
        'ID', 'Sample', 'Cohort',  # 'Subtype', 'ConfidenceLevel', 'Study'
    ]


def parse(file):

    cohort = os.path.basename(file).split('.')[0]

    samples = set()
    for m in TSVReader(file):
        _, sample, ref, alt, pos = m['#Uploaded_variation'].split('__')
        samples.add(sample)

    for i, sample in enumerate(samples):
        fields = [
            str(i),
            sample,
            cohort
        ]
        yield fields
