import csv
import gzip
import json
import os
from collections import defaultdict

import bgdata
import click
import numpy as np
from bgreference import refseq
from intervaltree import IntervalTree
from pyliftover import LiftOver

from intogen_core.parsers import DatasetError
from intogen_core.readers import TSVReader
from intogen_core.utils import out_open

CHROMOSOMES = set([str(c) for c in range(1, 23)] + ['X', 'Y'])

CHR_MAX = {}
for chr_ in CHROMOSOMES:
    CHR_MAX[chr_] = len(refseq('hg38', chr_, 1, -1))


class _NoLiftOver:

    def convert_coordinate(self, chr_, pos, strand):
        return [(chr_, pos)]


def liftover_factory(from_, to):
    if from_ == to:
        return _NoLiftOver()
    else:
        return LiftOver(from_, to_db=to, search_dir=os.environ['INTOGEN_DATASETS']+'/liftover')


def __none_to_string(value):
    if value is None:
        return "None"
    return value


def hypermutators_cutoff(snv_per_sample, cutoff):
    vals = list(snv_per_sample.values())
    iqr = np.subtract(*np.percentile(vals, [75, 25]))
    q3 = np.percentile(vals, 75)
    computed_cutoff = (q3 + 1.5 * iqr)
    cutoff = max(cutoff, computed_cutoff)
    return cutoff, computed_cutoff, set([k for k, v in snv_per_sample.items() if v > cutoff])


def filter_(file, genome, cutoff, stats):

    # Compute global stats
    donors, mut_per_sample, snv_per_sample, indel_per_sample = \
        defaultdict(set), defaultdict(int), defaultdict(int), defaultdict(int)
    for m in TSVReader(file):
        s = m['SAMPLE']
        if m['ALT_TYPE'] == 'snv':
            snv_per_sample[s] += 1
        elif m['ALT_TYPE'] == 'indel':
            indel_per_sample[s] += 1

        mut_per_sample[s] += 1
        donors[m['DONOR']].add(s)

    if len(snv_per_sample) < 1:
        raise DatasetError('No samples with SNVs')

    if len(donors) == 1 and ("None" in donors or "" in donors):
        # assume donor was not provided to bgparsers before
        donors = {}

    for d in donors:
        donors[d] = list(sorted(donors[d]))

    stats['donors'] = {__none_to_string(d): s for d, s in donors.items()}

    if None in donors:
        stats['warning_no_donor_id'] = "There is no donor ID"
    donors_with_multiple_samples = [d for d, s in donors.items() if len(s) > 1]
    if len(donors_with_multiple_samples) > 0:
        stats['warning_multiple_samples_per_donor'] = "{}".format(donors_with_multiple_samples)

    # We only want to use one sample per donor
    # this dictionary contains the samples that we'll skip
    multiple_donor_samples = []
    for d in donors_with_multiple_samples:
        multiple_donor_samples += donors[d][1:]
    multiple_donor_samples = set(multiple_donor_samples)

    cutoff, theorical_cutoff, hypermutators = hypermutators_cutoff(snv_per_sample, cutoff)
    stats['hypermutators'] = {
        'cutoff': cutoff,
        'computed_cutoff': theorical_cutoff,
        'hypermutators': list(hypermutators)
    }

    # Load coverage regions tree
    # TODO do it only on hg38
    coverage_tree = defaultdict(IntervalTree)
    with gzip.open(bgdata.get(f'intogen/coverage/{genome}'), 'rt') as fd:
        reader = csv.reader(fd, delimiter='\t')
        for i, r in enumerate(reader, start=1):
            coverage_tree[r[0]][int(r[1]):(int(r[2]) + 1)] = i

    # Load Somatic PONs file
    somatic_pon = set()
    somatic_pon_file = os.path.join(os.environ['INTOGEN_DATASETS'], 'others', 'somatic_pon_count_filtered.tsv.gz')
    with gzip.open(somatic_pon_file, 'rt') as fd:
        for line in fd:
            line = tuple(line.strip().split('\t'))
            somatic_pon.add(line)

    lo = liftover_factory(genome, 'hg38')
    # TODO add to docs that only hg19 and 38 are supported

    skipped = defaultdict(int)
    skip_chromosome, skip_chromosome_names = 0, set()
    skip_coverage, skip_coverage_positions = 0, []
    variants_by_sample = defaultdict(set)
    signature = defaultdict(int)

    counts = defaultdict(int)
    for v in TSVReader(file):
        counts['before'] += 1

        v['POSITION'] = int(v['POSITION'])

        # Force the strand to be positive
        v['STRAND'] = '+'

        if v['REF'] == v['ALT']:
            skipped['same_alt'] += 1
            continue

        # Skip hypermutators and multiple samples per donor
        if v['SAMPLE'] in hypermutators:
            skipped['hypermutators'] += 1
            continue

        # Skip multiple samples per donor
        if v['SAMPLE'] in multiple_donor_samples:
            skipped['multiple_samples_per_donor'] += 1
            continue

        if v['CHROMOSOME'] not in CHROMOSOMES:
            skip_chromosome_names.add(v['CHROMOSOME'])
            skip_chromosome += 1
            continue

        if v['CHROMOSOME'] in coverage_tree:
            if len(coverage_tree[v['CHROMOSOME']][v['POSITION']]) == 0:
                skip_coverage += 1
                skip_coverage_positions.append((v['SAMPLE'], v['CHROMOSOME'], v['POSITION']))
                continue

        # Check duplicates
        var_value = "{}:{}:{}>{}".format(v['CHROMOSOME'], v['POSITION'], v['REF'], v['ALT'])
        if var_value in variants_by_sample[v['SAMPLE']]:
            skipped['duplicated'] += 1
            continue
        variants_by_sample[v['SAMPLE']].add(var_value)

        # Skip variants that has an N
        if 'N' in v['ALT'] or 'N' in v['REF']:
            skipped['n_sequence'] += 1
            continue

        if v['ALT_TYPE'] == 'snv':
            # Compute signature and count mismatch
            ref = refseq(genome, v['CHROMOSOME'], v['POSITION'], size=1).upper()
            if ref != v['REF']:
                skipped['mismatch'] += 1
                continue

        # Liftover to hg38
        strand = '+' if 'STRAND' not in v else v['STRAND']

        # TODO rewrite this part to make it more meaningful
        lo_pos = lo.convert_coordinate(f"chr{v['CHROMOSOME']}", v['POSITION'] - 1, strand)

        if lo_pos is None or len(lo_pos) != 1 or lo_pos[0][0] is None:
            skipped['noliftover'] += 1
            continue

        lo_chr = lo_pos[0][0].replace('chr', '')
        if lo_chr not in CHROMOSOMES:
            skip_chromosome_names.add(lo_chr)
            skip_chromosome += 1
            continue

        lo_pos = lo_pos[0][1] + 1
        # check the position in in the chromosome
        # TODO check, because the old version compares against the other genome
        if lo_pos < 1 or lo_pos > CHR_MAX[lo_chr]:
            skipped['noliftover'] += 1
            continue

        v['CHROMOSOME'] = lo_chr
        v['POSITION'] = lo_pos

        # Skip variants that are in the somatic_pon_count_filtered.tsv.gz file
        var_value = (v['CHROMOSOME'], str(v['POSITION']), v['REF'], v['ALT'])
        if var_value in somatic_pon:
            skipped['somatic_pon'] += 1
            continue

        if v['ALT_TYPE'] == 'snv':
            # Compute signature and count mismatch
            ref = refseq('hg38', v['CHROMOSOME'], v['POSITION'] - 1, size=3).upper()
            alt = ''.join([ref[0], v['ALT'], ref[2]])
            if ref[1] != v['REF']:
                skipped['liftover_mismatch'] += 1
                continue
            signature_key = "{}>{}".format(ref, alt)
            signature[signature_key] = signature.get(signature_key, 0) + 1

            counts['snv'] += 1
        elif v['ALT_TYPE'] == 'indel':
            counts['indel'] += 1
        counts['after'] += 1

        yield v

    skipped['invalid_chromosome'] = (skip_chromosome, list(skip_chromosome_names))
    skipped['coverage'] = (skip_coverage, skip_coverage_positions)
    stats['skip'] = skipped
    stats['count'] = counts
    stats['signature'] = signature

    signature_count = sum(signature.values())
    if signature_count > 0:
        stats['probabilities'] = {k: v / signature_count for k, v in signature.items()}
    else:
        stats['probabilities'] = signature

    mismatches = stats['skip']['mismatch']
    snvs = stats['count']['snv']
    ratio_mismatch = (mismatches / snvs) if snvs > 0 else 0
    if ratio_mismatch > 0.1:
        raise DatasetError(f'There are {mismatches} of {snvs} genome reference mismatches. More than 10%.')
    elif ratio_mismatch > 0.05:
        stats["warning_genome_reference_mismatch"] = f"There are {mismatches} of {snvs} genome reference mismatches."
    same_alt = stats['skip']['same_alt']
    if same_alt > 0:
        stats["warning_same_alternate"] = f"There are {same_alt} entries with same reference and alternate"

    lo_mismatches = stats['skip']['liftover_mismatch']
    ratio_lo_mismatch = (lo_mismatches / snvs) if snvs > 0 else 0
    if ratio_lo_mismatch > 0.3:
        raise DatasetError(f'There are {lo_mismatches} of {snvs} genome reference mismatches after liftover. More than 30%.')

    count_after = stats['count']['after']
    if count_after == 0:
        raise DatasetError('There are no variants after filtering')

    duplicated = stats['skip']['duplicated']
    if duplicated > 0:
        stats["warning_duplicated_variants"] = f"There are {duplicated} duplicated variants"

    ns = stats['skip']['n_sequence']
    if ns > 0:
        stats["warning_n_sequence"] = f"There are {ns} variants with a 'N' in the reference or alternate sequence"


def process(infile, outfile, genome, cutoff):
    stats = {}

    try:

        with out_open(outfile, 'wt') as fd:
            writer = None
            for row in filter_(infile, genome, cutoff, stats):
                if writer is None:
                    fields = [f for f in row.keys() if f not in ['DATASET', 'PLATFORM', 'GENOMEREF']]
                    writer = csv.DictWriter(fd, fieldnames=fields, delimiter='\t', extrasaction='ignore')
                    writer.writeheader()
                writer.writerow(row)

    except DatasetError as e:
        # remove the file
        os.unlink(outfile)
        raise
    else:
        # TODO check if stats are useful, or we should just print
        # write stats to file
        with open(outfile + '.stats.json', 'w') as fd:
            json.dump(stats, fd, indent=4)


@click.command()
@click.option('-i', '--input', type=click.Path(exists=True), required=True)
@click.option('-o', '--output', type=click.Path(), required=True)
@click.option('-g', '--genome', type=click.Choice(['hg19', 'hg38']), default='hg38')
@click.option('-c', '--cutoff', default=1000)
def cli(input, output, genome, cutoff):
    process(input, output, genome, cutoff)


if __name__ == '__main__':
    cli()
