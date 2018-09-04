import os
import csv
import logging

from .base import Filter
from tasks.base import valid_consequence

logger = logging.getLogger(__name__)


class VepFilter(Filter):

    KEY = "vep"

    def __init__(self, parent):
        super().__init__(parent)

        with open(os.path.join(os.environ['INTOGEN_DATASETS'], 'selected_ensembl_proteins.tsv')) as fd:
            self.proteins = set([r[2] for r in csv.reader(fd, delimiter='\t')])

        with open(os.path.join(os.environ['INTOGEN_DATASETS'], 'selected_ensembl_proteins.tsv')) as fd:
            self.genes = set([r[1] for r in csv.reader(fd, delimiter='\t')])

    def run(self, group_key, group_data):

        # To store errors and statistics
        self.stats[group_key] = {}

        self.KEY = "vep{}{}".format(os.path.sep, group_key)

        genes = {}
        consequence = {}
        chromosomes = {}

        count_before = 0
        count_skip_consequence = 0
        count_after = 0

        skip_genes = set()
        process_genes = set()


        for v in self.parent.run(group_key, group_data):
            count_before += 1

            if not valid_consequence(v['Consequence']):
                count_skip_consequence += 1
                continue

            if v['ENSP'] not in self.proteins:
                if v['Gene'] in self.genes:
                    # Selected gene without matching protein id
                    skip_genes.add(v['Gene'])

                continue

            process_genes.add(v['Gene'])

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
            'before': count_before,
            'skip_consequence': count_skip_consequence
        }

        self.stats[group_key]['ratio_missense'] = consequence.get('missense_variant', 0) / consequence.get('synonymous_variant', 0) if 'synonymous_variant' in consequence else None

        orphan_genes = skip_genes - process_genes
        if len(orphan_genes) > 0:
            self.stats[group_key]['orphan_genes'] = [g for g in orphan_genes]
            self.stats[group_key]["warning_orphan_genes"] = "There are {} orphan genes at {}".format(len(orphan_genes), group_key)

        if count_after == 0:
            self.stats[group_key]["error_no_entries"] = "There is no VEP output"

        if len(chromosomes) < 14:
            self.stats[group_key]["error_few_chromosomes_with_mutations"] = "There are only {} chromosomes with mutations at {}".format(len(chromosomes), group_key)
        elif len(chromosomes) < 23:
            self.stats[group_key]["warning_few_chromosomes_with_mutations"] = "There are only {} chromosomes with mutations at {}".format(len(chromosomes), group_key)

