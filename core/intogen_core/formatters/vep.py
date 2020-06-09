
from intogen_core.readers import TSVReader

HEADER = None


def parse(file):
    i = 0
    for m in TSVReader(file):
        # TODO what happens with positions of indels? VEP requires other format
        fields = [
            m['CHROMOSOME'],
            m['POSITION_HG38'],
            m['POSITION_HG38'],
            f"{m['REF']}/{m['ALT']}",
            m['STRAND'],
            f"I{i:010d}__{m['SAMPLE']}__{m['REF']}__{m['ALT']}"
        ]
        yield fields
        i += 1

# TODO the original dataset includes an extra line as "{}__{}__{}__{}".format(identifier, mut['SAMPLE'], mut['REF'], mut['ALT']). Is needed?