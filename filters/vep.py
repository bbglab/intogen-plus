import os
import logging

from .base import Filter

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
        chromosomes = {}

        count_before = 0
        count_after = 0

        for v in self.parent.run(group_key, group_data):
            count_before += 1

            if v['CANONICAL'] != "YES":
                continue

            # Remove multiple consequences
            v['Consequence'] = v['Consequence'].split(',')[0]

            chromosome = v['Location'].split(":")[0]
            chromosomes[chromosome] = chromosomes.get(chromosome, 0) + 1
            count_after += 1
            consequence[v['Consequence']] = consequence.get(v['Consequence'], 0) + 1
            genes[v['SYMBOL']] = genes.get(v['SYMBOL'], 0) + 1

            yield v

        self.stats[group_key]['consequence'] = consequence
        self.stats[group_key]['chromosomes'] = chromosomes
        self.stats[group_key]['genes'] = genes
        self.stats[group_key]['count'] = {
            'after': count_after,
            'before': count_before
        }
        self.stats[group_key]['ratio_missense'] = consequence.get('missense_variant', 0) / consequence.get('synonymous_variant', 0) if 'synonymous_variant' in consequence else None

        if count_after == 0:
            self.stats[group_key]["error_no_entries"] = "There is no VEP output"

        if len(chromosomes) < 14:
            self.stats[group_key]["error_few_chromosomes_with_mutations"] = "There are only {} chromosomes with mutations".format(len(chromosomes))
        elif len(chromosomes) < 23:
            self.stats[group_key]["warning_few_chromosomes_with_mutations"] = "There are only {} chromosomes with mutations".format(len(chromosomes))

