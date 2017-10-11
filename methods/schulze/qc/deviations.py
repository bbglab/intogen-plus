# Import modules
from collections import namedtuple

import os
import pandas as pd
import numpy as np
import logging
from scipy.stats import linregress

from .drivers import CGC_GENES_PER_TUMOR

# Configure the logger
logger = logging.getLogger(__name__)

Columns = namedtuple('Parser', 'GENE_ID SYMBOL PVALUE QVALUE')


class Deviation:
    def __init__(self, df=None, description=''):
        self.df = df
        self.description = description
        self.parser = Columns(GENE_ID='GENE_ID', SYMBOL='SYMBOL', PVALUE='PVALUE', QVALUE='QVALUE')

    def deviation_from_null(self, half):
        """Calculate the deviation from the null hypothesis
        :param half: boolean, if True analyze only half of the distribution of pvalues
        :return: float, deviation from the null
        """
        # Calculate the minimum p-value
        min_pvalue = float(np.min(self.df[self.df[self.parser.PVALUE] > 0][[self.parser.PVALUE]]))

        obs_pvalues = sorted(
            self.df[self.parser.PVALUE].map(lambda x: -np.log10(x) if x > 0 else -np.log10(min_pvalue))
        )
        exp_pvalues = -np.log10(np.arange(1, len(self.df) + 1) / float(len(self.df)))
        exp_pvalues.sort()

        # Get half distribution
        if half:
            obs_pvalues = obs_pvalues[:len(self.df) // 2]
            exp_pvalues = exp_pvalues[:len(self.df) // 2]

        if len(obs_pvalues) < 10:
            logger.warning('Cannot calculate the deviation in {} (< 10 points)'.format(self.description))
            return {'deviation': np.nan, 'slope': np.nan}

        # Calculate the average square difference
        deviation = np.mean([(i - j) ** 2 for i, j in zip(obs_pvalues, exp_pvalues)])

        # Calculate the slope
        try:
            linear_regression = linregress(exp_pvalues, obs_pvalues)
        except ValueError as err:
            return {'deviation': np.nan, 'slope': np.nan}
        return {'deviation': deviation, 'slope': linear_regression.slope}

    @staticmethod
    def get_weight(i, weight):
        """Calculate the weight of an ith position
        :param i: the ith ranking
        :param weight: the type of weighting [log, normal]. Default normal.
        :return: the weight of the ith position
        """
        if weight == "log":
            return 1.0 / np.log2(i + 2)
        if weight == "normal":
            return 1.0 / i

    @staticmethod
    def calculate_percentage_cgc(ranking):
        """Calculate the percentage of CGC genes in the input list
        :param ranking: the input list of the ranked genes
        :return: percentage of cgc genes in the list
        """
        n = float(sum([1.0 if gene in CGC_GENES_PER_TUMOR['PANCANCER'] else 0.0 for gene in ranking]))
        return n / len(ranking)

    def get_maximum_area(self, position, weight):
        """Returns the maxmimum theoric weighted area for that position
        :param position: the position of the max value
        :param weight: type of normalization
        :return: maximum theoric area
        """
        return sum([1 * self.get_weight(i, weight) for i in range(position + 1)])

    def calculate_areas(self, up_to=100):
        """Calculate the ranking absolute and relative areas. The relative area is
        the area under the curve of the CGC enrichment of a given rankings normalized
        by the maximum reachable area by the number of ranked genes.
        :param up_to: int, maximum number of genes to consider
        :return: dictionary representing the two areas, absolute and relative
        """
        # positive = self.df[self.df[self.parser.QVALUE] < 0.1]
        # up_to = min(up_to, len(positive))
        self.df.sort_values(by=self.parser.PVALUE, ascending=True, inplace=True)
        ranking = self.df[:up_to][self.parser.SYMBOL].tolist()
        xticks = range(len(ranking))
        area = 0.0
        for i in xticks:
            weight_i = self.get_weight(i, weight='log')
            x_i = self.calculate_percentage_cgc(ranking[0: i + 1])
            area += x_i * weight_i
        max_area = self.get_maximum_area(len(ranking), weight='log')
        return {'absolute': area, 'relative': area / max_area}

