
import csv

import click

from intogen_core.formatters import fml, clustl, vep, signature, dndscv, smregions, cbase,  mutpanning_muts, mutpanning_samples, deconstructsigs, hotmaps
from intogen_core.utils import out_open

MAP = {
    'fml': fml,
    'clustl': clustl,
    'vep': vep,
    'signature': signature,
    'dndscv': dndscv,
    'smregions': smregions,
    'cbase': cbase,
    'mutpanning-mutations': mutpanning_muts,
    'mutpanning-samples': mutpanning_samples,
    'deconstructsigs': deconstructsigs,
    'hotmaps': hotmaps
}


@click.command()
@click.option('-i', '--input', type=click.Path(exists=True), required=True)
@click.option('-o', '--output', type=click.Path(), required=True)
@click.option('-f', '--format', type=click.Choice(list(MAP.keys())))
def cli(input, output, format):

    module = MAP[format.lower()]

    with out_open(output, 'wt') as fd:
        writer = csv.writer(fd, delimiter='\t')
        if module.HEADER is not None:
            writer.writerow(module.HEADER)
        for row in module.parse(input):
            writer.writerow(row)


if __name__ == '__main__':
    cli()
