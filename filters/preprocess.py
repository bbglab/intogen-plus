import os
import csv
import gzip
import logging
import numpy as np

from .base import Filter
from collections import Counter, defaultdict
from intervaltree import IntervalTree

logger = logging.getLogger(__name__)


class PreprocessFilter(Filter):

    # Minimum cutoff
    MIN_CUTOFF = 1000
    CHROMOSOMES = set(list(range(1, 23)) + ['X', 'Y'])

    def __init__(self, source):
        self.source = source

    def run(self, group_key, group_data):

        # Find Hypermutators Samples
        sample_muts = Counter(
            [m['SAMPLE'] for m in self.source.run(group_key, group_data) if m['ALT_TYPE']=='snp']
        )

        vals = list(sample_muts.values())
        if len(vals) == 0:
            return

        iqr = np.subtract(*np.percentile(vals, [75, 25]))
        q3 = np.percentile(vals, 75)
        cutoff = max(self.MIN_CUTOFF, (q3 + 1.5 * iqr))
        hypermutators = set([k for k, v in sample_muts.items() if v > cutoff])
        if len(hypermutators) > 0:
            logger.info("[QC] {} HYPERMUTATORS at {}:  {}".format(group_key, cutoff, ", ".join(["{} = {}".format(h, sample_muts[h]) for h in hypermutators])))

        # Load coverage regions tree
        regions_file = os.environ['COVERAGE_REGIONS']
        coverage_tree = defaultdict(IntervalTree)
        with gzip.open(regions_file, 'rt') as fd:
            reader = csv.reader(fd, delimiter='\t')
            for i, r in enumerate(reader, start=1):
                coverage_tree[r[0]][int(r[1]):(int(r[2]) + 1)] = i

        # Read variants
        for v in self.source.run(group_key, group_data):

            # Skip hypermutators
            if v['SAMPLE'] in hypermutators:
                continue

            if v['CHROMOSOME'] not in self.CHROMOSOMES:
                continue

            if v['CHROMOSOME'] in coverage_tree:
                if len(coverage_tree[v['CHROMOSOME']][v['POSITION']]) == 0:
                    logger.info("[QC] {} LOW COVERAGE: {} at {}:{}".format(group_key, v['SAMPLE'], v['CHROMOSOME'], v['POSITION']))
                    continue
            yield v

