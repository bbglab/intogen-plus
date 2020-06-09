import os
import random

import numpy as np
import pandas as pd



class Evaluation_Enrichment:

    @staticmethod
    def randomize_subset(cgc, percentage=1.0):
        total = len(cgc)
        selected = int(total * percentage)
        return set(random.sample(list(cgc), selected))

    @staticmethod
    def load_cgc(percentage=1.0):
        cgc = pd.read_csv(os.path.join(
            os.environ["INTOGEN_DATASETS"], 'cgc', "cancer_gene_census_parsed.tsv"
        ), sep="\t")
        cgc = set(cgc["Gene Symbol"].values)
        if percentage < 1.0:
            cgc_f = Evaluation_Enrichment.randomize_subset(cgc, percentage)
        else:
            cgc_f = cgc
        return cgc_f

    def __init__(self, percentage):

        self.type_evaluation = "CGC_Enrichment"
        self.cgc = Evaluation_Enrichment.load_cgc(percentage)

    def get_weight(self, i, weight):
        """
        :param i: the ith ranking
        :param weight: the type of weighting [log,normal]. Default normal.
        :return: the weight of the ith position
        """
        if weight == "log":
            return 1.0 / np.log2(i+2)
        if weight == "normal":
            return 1.0 / i

    def calculate_percentage_cgc(self, ranking):
        """
        Calculate the percentage of CGC genes in the input list
        :param ranking: the input list of the ranked genes

        :return: percentage of cgc genes in the list
        """
        n = float(sum([1.0 if gene in self.cgc else 0.0 for gene in ranking]))
        return n / len(ranking)

    def evaluate_enrichment_method(self, ranking, weight="log", ranking_limit=40):
        """
        It calculates the area under the curve of the CGC enrichment of the given rankings.
        :param ranking: sorted list of genes
        :param weight: normalization of the weight. [log,normal] default: log
        :param cancer_drivers: gold-standardt of cgc drivers. Used to calculate the enrichement.
        :param ranking_limit: limit to calculate the area under the curve. Default: top-40
        :return The weighted area under the curve of the CGC enrichment
        """
        if len(ranking) > ranking_limit:
            ranking = ranking[0:ranking_limit]
        xticks = range(len(ranking))

        area = 0.0
        for i in xticks:

            weight_i = self.get_weight(i, weight)

            x_i = self.calculate_percentage_cgc(ranking[0:i+1])

            area = area + x_i*weight_i
        return area

    def get_maximum_area(self, position, weight):
        """
        Returns the maxmimum theoric weighted area for that position
        :param position: the position of the max value
        :param weight: type of normalization
        :return: maximum theoric area
        """
        maxv = sum([1*self.get_weight(i, weight) for i in range(position+1)])
        return maxv

    def evaluate_enrichment_method_relative(self, ranking, weight="log"):
        """
        It calculates the RELATIVE area under the curve of the CGC enrichment of the given rankings.
        It normalizes by the maximum reachable area by the number of ranked genes.
        :param ranking: sorted list of genes
        :param weight: normalization of the weight. [log,normal] default: log
        :param cancer_drivers: gold-standardt of cgc drivers. Used to calculate the enrichement.
        :param ranking_limit: limit to calculate the area under the curve. Default: top-40
        :return The weighted area under the curve of the CGC enrichment
        """

        xticks = range(len(ranking))

        area = 0.0
        for i in xticks:

            weight_i = self.get_weight(i, weight)

            x_i = self.calculate_percentage_cgc(ranking[0:i+1])

            area = area + x_i*weight_i
        max_area = self.get_maximum_area(len(ranking), weight)

        return float(area)/float(max_area)  # Max value 1

    def get_list_ranking(self, d_results):
        """
        Given a dictionary with the ranking of the genes returns a list with keeping the ranking of the genes.
        :param d_results: dictionary of format {gene:ranking,gene2:ranking2,gene3:ranking3....}
        :return: a list with the sorted ranking of genes
        """

        return [x[0] for x in sorted(d_results.items(), key=lambda x: (x[1], x[0]))]

    def calculate_area_cancer(self, d_results, type_method):
        """
        :param d_results: dictionary of results including all cancer types and methods
        :return dictionary with the area per cancer type per method
        """
        d_area = {}

        for cancer in d_results:
            d_area[cancer] = {}
            for method in d_results[cancer]:
                list_ranking = self.get_list_ranking(d_results[cancer][method])
                if type_method == "absolute":
                    d_area[cancer][method] = self.evaluate_enrichment_method(list_ranking)
                else:
                    d_area[cancer][method] = self.evaluate_enrichment_method_relative(list_ranking)

        return d_area

    def calculate_area(self, d_results, type_method="absolute"):
        """
        :param  d_results: dictionary of results for a given cancer type
        :return  dictionary with the area per method
        """
        d_area = {}
        for method in d_results.keys():

            list_ranking = self.get_list_ranking(d_results[method])
            if type_method == "absolute":
                d_area[method] = self.evaluate_enrichment_method(list_ranking)
            else:
                d_area[method] = self.evaluate_enrichment_method_relative(list_ranking)

        return d_area

    def calculate_area_list(self, list_ranking, type_method="absolute"):
        """
        :param list_ranking: list with the genes sorted by ranking
        :return the value of the AUC for this ranking
        """
        if type_method == "absolute":
            return self.evaluate_enrichment_method(list_ranking)
        else:
            return self.evaluate_enrichment_method_relative(list_ranking)

    @staticmethod
    def get_report_area(dict_results):
        """
        :param dict_results: the input dictionary of weighted AUC
        :return: a DataFrame with the information
        """
        l_results = []
        for cancer in dict_results.keys():
            for method in dict_results[cancer]:
                l_results.append([cancer, method, dict_results[cancer][method]])

        return pd.DataFrame(l_results, columns=["Cancer_Type", "Method", "Weighted_AUC"])
