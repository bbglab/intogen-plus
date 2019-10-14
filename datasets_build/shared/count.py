import csv
import gzip
import json
from collections import defaultdict
from itertools import product
import pickle

import bglogs
import click
import numpy as np
from bgreference import refseq
from bgvep.readers import BGPack


def open_(file):
    if file.endswith('.gz'):
        return gzip.open(file, 'rt')
    else:
        return open(file)


def triplets(sequence):
    iterator = iter(sequence)

    n1 = next(iterator)
    n2 = next(iterator)

    for n3 in iterator:
        yield n1 + n2 + n3
        n1 = n2
        n2 = n3


def compute(regions, genome, vep):

    gene_triplets = defaultdict(lambda: defaultdict(int))
    gene_conseq = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    counter = 0

    with BGPack(genome, vep) as conseq_reader:
        for gene, chr_, start, stop in regions:

            counter += 1
            if counter % 10000 == 0:
                bglogs.info('---%d' % counter)

            bglogs.debug('Region %s, %s, %d, %d' % (gene, chr_, start, stop))

            try:
                seq = refseq(genome, chr_, start-1, stop-start+3)
                bglogs.debug('Sequence %s' % seq)
            except (ValueError, RuntimeError):
                bglogs.warning('Error in ref for CHR: %s positions: %d:%d' % (chr_, start, stop))
                continue
            for triplet, v in zip(triplets(seq), conseq_reader.get(chr_, start, stop)):
                gene_triplets[gene][triplet] += 1
                ref = triplet[1]
                pos, consequences = v
                bglogs.debug('Triplet %s at %s has consequences %s' % (triplet, pos, consequences))
                for i, alt in enumerate('ACGT'):
                    if alt == ref:
                        continue
                    gene_conseq[gene][triplet+'>'+alt][consequences[i]] += 1

    return gene_triplets, gene_conseq


def save_counts(counts, file):
    with gzip.open(file, 'wt') as fd:
        json.dump(counts, fd, indent=4)
    bglogs.info('Saved to {}'.format(file))


def pickle_counts(counts, file):
    with gzip.open(file, 'wb') as fd:
        pickle.dump(counts, fd)
    bglogs.info('Saved to {}'.format(file))


def load_regions(file):
    with open_(file) as fd:
        for row in csv.DictReader(fd, delimiter='\t',
                                  fieldnames=['CHR', 'START', 'STOP', 'STRAND', 'GENE', 'TRANSCRIPT', 'SYMBOL']):
            gene = row['SYMBOL']
            chr_ = row['CHR']
            start = int(row['START'])
            stop = int(row['STOP'])

            yield gene, chr_, start, stop


def context_generator():
    subs = list(filter(lambda x: x[0] != x[1], map(lambda x: ''.join(x), product(list('ACGT'), repeat=2))))
    for s in sorted(subs):
        for f in sorted(product(list('ACGT'), repeat=2)):
            yield f[0] + s[0] + f[1] + '>' + s[1]


def compress(d):
    genes = sorted(d.keys())
    contexts = list(context_generator())
    csqn_types = ['synonymous_variant', 'missense_variant', 'stop_gained_variant',
                  'splice_acceptor_variant', 'splice_donor_variant']
    count_matrix = np.zeros((len(genes), len(contexts), len(csqn_types)))

    for i, gene in enumerate(genes):
        for j, context in enumerate(contexts):
            csqn_count = d[gene].get(context, {})
            for k, csqn_type in enumerate(csqn_types):
                count_matrix[i, j, k] = csqn_count.get(csqn_type, 0)

    compressed = {'genes': genes, 'contexts': contexts, 'csqn_types': csqn_types, 'matrix': count_matrix}

    return compressed


def run(regions_file, triplets_file, conseq_file, genome, vep):
    gene_triplets, gene_consq = compute(load_regions(regions_file), genome, vep)
    save_counts(gene_triplets, triplets_file)
    gene_consq = compress(gene_consq)
    pickle_counts(gene_consq, conseq_file)


@click.command()
@click.option('-r', '--regions', 'regions_file', help='Path to regions file', type=click.Path(), required=True)
@click.option('-t', '--triplets', 'triplets_file', help='Output path for triplets', type=click.Path(), required=True)
@click.option('-c', '--consequence', 'conseq_file', help='Output path for consequences', type=click.Path(), required=True)
@click.option('-v', '--vep', 'vep_version', help='VEP version', required=True)
@click.option('-g', '--genome', 'genome_build', help='Reference genome', required=True)
@click.option('--debug', is_flag=True)
def cli(regions_file, triplets_file, conseq_file, vep_version, genome_build, debug):
    bglogs.configure(debug=debug)

    run(regions_file, triplets_file, conseq_file, genome_build, vep_version)


if __name__ == '__main__':
    cli()
