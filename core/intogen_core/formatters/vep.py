
from intogen_core.readers import TSVReader

HEADER = None


def parse(file):
    i = 0
    for m in TSVReader(file):

        start = end = int(m["POSITION"])
        if m["REF"] == "-":  # is an insertion
            start = end + 1
        elif m["ALT"] == "-":  # is a deletion
            end = start + len(m["REF"]) - 1
        else:  # snv/mnv
            end = start + len(m["ALT"]) - 1

        fields = [
            m['CHROMOSOME'],
            f"{start}",
            f"{end}",
            f"{m['REF']}/{m['ALT']}",
            m['STRAND'],
            f"I{i:010d}__{m['SAMPLE']}__{m['REF']}__{m['ALT']}__{m['POSITION']}"
        ]
        yield fields
        i += 1
