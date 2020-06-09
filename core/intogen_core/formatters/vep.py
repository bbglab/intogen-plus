
from intogen_core.readers import TSVReader

HEADER = None


def parse(file):
    i = 0
    for m in TSVReader(file):
        # TODO what happens with positions of indels? VEP requires other format
        fields = [
            m['CHROMOSOME'],
            m['POSITION'],
            m['POSITION'],
            f"{m['REF']}/{m['ALT']}",
            m['STRAND'],
            f"I{i:010d}__{m['SAMPLE']}__{m['REF']}__{m['ALT']}"
        ]
        yield fields
        i += 1
