import os
import csv
import gzip
import logging
import numpy as np

from bgreference import hg19
from .base import Filter
from collections import Counter, defaultdict
from intervaltree import IntervalTree

logger = logging.getLogger(__name__)


class VepFilter(Filter):

    KEY = "vep"

    def __init__(self, parent):
        super().__init__(parent)

    def run(self, group_key, group_data):

        # To store errors and statistics
        self.stats[group_key] = {}

        self.KEY = "vep{}{}".format(os.path.sep, group_key)

        genes = {}
        consequence = {}

        count_before = 0
        count_after = 0

        for v in self.parent.run(group_key, group_data):
            count_before += 1

            if v['CANONICAL'] != "YES":
                continue

            count_after += 1
            consequence[v['Consequence']] = consequence.get(v['Consequence'], 0) + 1
            genes[v['SYMBOL']] = genes.get(v['SYMBOL'], 0) + 1

            yield v

        self.stats[group_key]['consequence'] = consequence
        self.stats[group_key]['genes'] = genes
        self.stats[group_key]['count'] = {
            'after': count_after,
            'before': count_before
        }

        if count_after == 0:
            self.stats[group_key]["error_no_entries"] = "There is no VEP output"
