import csv

import click
from bgparsers import readers
from bgreference import hg38


def generate_all_possible_snvs(regions):

    for region in readers.elements(regions):
        chr_, start, end = region['CHROMOSOME'], region['START'], region['END']
        for i, ref in enumerate(hg38(chr_, start, end - start+1), start=start):
            if ref == 'N':
                continue
            for alt in 'ACGT':
                if ref == alt:
                    continue
                yield chr_, i, ref, alt



@click.command()
@click.option('-r', '--regions', type=click.Path(exists=True), required=True)
@click.option('-o', '--output', type=click.Path(), required=True)
def cli(regions, output):

    with open(output, 'wt') as fd:
        writer = csv.writer(fd, delimiter='\t')
        for chr_, pos, ref, alt in generate_all_possible_snvs(regions):
            writer.writerow(
                [
                    chr_,
                    f"{pos}",
                    f"{pos}",
                    f"{ref}/{alt}",
                ]
            )


if __name__ == '__main__':
    cli()