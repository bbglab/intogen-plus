import csv
import json
import os

import click

from intogen_core.parsers import DatasetError
from intogen_core.readers import TSVReader
from intogen_core.utils import out_open


def filter_(file, stats):

    synonymous = 0

    for v in TSVReader(file):
        if v['Consequence'] == 'synonymous_variant':
            synonymous += 1
        else:
            yield v
    if synonymous > 0:
        stats['warning_synonymous_mutations'] = f"There are {synonymous} synonymous mutations"


def process(infile, outfile):
    stats = {}

    try:
        with out_open(outfile, 'wt') as fd:
            writer = None
            for row in filter_(infile, stats):
                if writer is None:
                    fields = [f for f in row.keys()]
                    writer = csv.DictWriter(fd, fieldnames=fields, delimiter='\t')
                    writer.writeheader()
                writer.writerow(row)

    except DatasetError as e:
        # remove the file
        os.unlink(outfile)
        raise
    else:
        # write stats to file
        with open(outfile + '.stats.json', 'w') as fd:
            json.dump(stats, fd, indent=4)


@click.command()
@click.option('-i', '--input', type=click.Path(exists=True), required=True)
@click.option('-o', '--output', type=click.Path(), required=True)
def cli(input, output):
    process(input, output)


if __name__ == '__main__':
    cli()
