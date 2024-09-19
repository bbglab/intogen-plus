import csv
import json
import os
from collections import defaultdict

import click

from intogen_core.parsers import DatasetError
from intogen_core.readers import TSVReader
from intogen_core.utils import out_open


VALID_CONSEQUENCES = {
        "transcript_ablation", "splice_donor_variant", 
        "splice_acceptor_variant", "stop_gained",
        "frameshift_variant", "stop_lost", "start_lost",
        "initiator_codon_variant", "transcript_amplification",
        "feature_elongation", "feature_truncation",
        "inframe_insertion", "inframe_deletion",
        "protein_altering_variant",
        "missense_variant",
        "splice_donor_5th_base_variant",
        "splice_region_variant", "incomplete_terminal_codon_variant",
        "start_retained_variant", "stop_retained_variant",
        "synonymous_variant", "coding_sequence_variant"
    }



def valid_consequence(consequences):
    csq_set = {c for c in consequences.split(',')}
    return len(csq_set & VALID_CONSEQUENCES) > 0


def filter_(file, stats):

    counts = defaultdict(int)
    chromosomes, genes, consequence = defaultdict(int), defaultdict(int), defaultdict(int)
    genes_skipped, genes_processed = set(), set()
    data = defaultdict(list)

    for v in TSVReader(file):
        counts['before'] += 1

        if not valid_consequence(v['Consequence']):
            counts['skip_consequence'] += 1
            continue

        if v['MANE_SELECT'] == '-':
            genes_skipped.add(v['Gene'])
            continue
        genes_processed.add(v['Gene'])

        # Remove multiple consequences
        # TODO shouldn't we get the first valid consequence?
        v['Consequence'] = v['Consequence'].split(',')[0]

        chromosome = v['Location'].split(":")[0]
        chromosomes[chromosome] += 1
        counts['after'] += 1
        consequence[v['Consequence']] += 1
        genes[v['SYMBOL']] += 1

        # Store all the lines matching a single position in a given patient
        # to get rid of positions that lie in more than one gene
        sample = v['#Uploaded_variation'].split('__')[1]
        data[(v['Location'], sample)].append(v)

    if counts['after'] == 0:
        raise DatasetError('There is no VEP output')

    synonymous_variants = consequence.get('synonymous_variant', 0)

    if synonymous_variants < 5:
        raise DatasetError(f'There are only {synonymous_variants} synonymous variants')

    if synonymous_variants == 0:
        ratio_missense = None
    else:
        ratio_missense = (consequence.get('missense_variant', 0) /
                          synonymous_variants)

    if ratio_missense is None:
        raise DatasetError("There are no synonymous variants")
    elif ratio_missense > 10:
        raise DatasetError("The ratio of missense / synonymous variants is too high")

    multiple_matches = {}
    for k, v in data.items():
        # If a (position , sample) (k) lie in more than one gene,
        # the first one in alphabetical order is selected and the
        # other are discarded
        if len(v) > 1:
            multiple_matches[k] = sorted([x['SYMBOL'] for x in v])
        yield sorted(v, key=lambda x: x['SYMBOL'])[0]

    stats['consequence'] = consequence
    stats['chromosomes'] = chromosomes
    stats['genes'] = genes
    stats['count'] = counts

    stats['ratio_missense'] = ratio_missense

    genes_orphan = genes_skipped - genes_processed
    if len(genes_orphan) > 0:
        stats['orphan_genes'] = list(genes_orphan)
        stats['warning_orphan_genes'] = f"There are {len(genes_orphan)} orphan genes"

    if len(chromosomes) < 23:
        stats['warning_few_chromosomes_with_mutations'] = f"There are only {len(chromosomes)} chromosomes with mutations"

    if len(multiple_matches) > 0:
        msg = ['{} in {} mapped to {}. Discarded {}'.format(k[0], k[1], v[0], ','.join(set(v[1:]))) for k, v in
               multiple_matches.items()]
        stats['warning_mutations_match_various_genes'] = msg


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
