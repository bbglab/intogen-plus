import os

import numpy as np
import pandas as pd

from intogen_combination.schulze_election import combination_ranking


def get_voters(gene, d):
    """
    Args:
        gene: str: gene symbol
        d: dict mapping method into dict mapping gene into rank
    Returns:
        methods: list of methods betting on symbol
        best_methods: highest bidding methods
        best_rank: best rank among voting methods
        ranks: list of ranks from those methods betting on symbol
    """

    d_rank = {}
    for method in d.keys():
        if gene in d[method]:
            d_rank[method] = d[method][gene]
        else:
            continue
    methods = list(d_rank.keys())
    if len(methods) == 0:
        print(gene)
    sorted_methods = sorted(d_rank.items(), key=lambda x: (x[1], x[0]))
    try:
        best_rank = sorted_methods[0][1]
        best_methods = [k for k, v in sorted_methods if v == best_rank]
        ranks = list(d_rank.values())
    except:
        best_rank = None
        best_methods = None
        ranks = None
    return methods, best_methods, best_rank, ranks


def output_to_dataframe(ranking, d_results):
    """
    Args:
        ranking: ranking dict: dict mapping candidates to ranks
        d_results: dict mapping methods to a ranking dict
    Returns:
        dataframe encoding summary information
    """

    l_info = []
    cgc = pd.read_csv(os.path.join(os.environ['INTOGEN_DATASETS'], "cgc", "cancer_gene_census_parsed.tsv"), sep="\t")
    cancer_drivers = cgc['Gene Symbol'].unique()

    for gene, rk in ranking.items():
        methods, best_methods, best_rank, ranks = get_voters(gene, d_results)
        try:
            median_rank = np.median(ranks)
            best_rank = min(ranks)
        except:
            median_rank = None
            best_rank = None
        l_info.append(
            [
                gene,
                rk,
                (gene in cancer_drivers),
                median_rank,
                best_rank,
                len(methods),
                ",".join(best_methods),
                ",".join(methods)
            ]
        )

    df_info = pd.DataFrame(l_info, columns=["SYMBOL", "RANKING", "CGC", "Median_Ranking", "Best_Ranking",
                                            "Total_Bidders", "Highest_Bidder", "All_Bidders"])
    return df_info


class Ballot(object):
    """
    dict mapping voter to dict mapping candidates to valid ranks
    ties are allowed
    """
    def __init__(self, ballot):
        self.dict = ballot

    def get_voter(self):
        return list(self.dict.keys()).pop()

    def get_candidates(self):
        voter = self.get_voter()
        return list(self.dict[voter].keys())

    def get_ranks(self):
        voter = self.get_voter()
        return list(self.dict[voter].values())

    def validate(self):
        a = None
        for i, rank in enumerate(sorted(self.get_ranks())):
            if a != rank:
                a = rank
                if rank != i + 1:
                    raise ValueError


def chunkizate(l, n_chunks):
    """
    Args:
        l: list
        n_chunks: int: number of chunks
    Returns:
        chunk_list = list of lists
    """
    n = len(l)
    q = n // n_chunks
    r = n % n_chunks
    chunk_list = []
    for i in range(n_chunks):
        if (r != 0) and (i < r):
            chunk_list.append(l[q*i: q*(i+1)] + [l[i-r]])
        else:
            chunk_list.append(l[q*i: q*(i+1)])
    return chunk_list


def strongest_paths_by_chunk(all_candidates, spath, chunk):
    """
    Args:
        chunk: list: list of candidates
        all_candidates: list: list of comprising all candidates
        spath: dict: dict mapping candidates to dict mapping candidates to strength of strongest path from
                     primary key to secondary key.
    Returns:
        updated spath: dict: dict mapping candidates to dict mapping candidates to strength of strongest path from
                     primary key to secondary key.
    """
    for i in chunk:
        for j in all_candidates:
            if i != j:
                for k in all_candidates:
                    if (i != k) and (j != k):
                        spath[i][j] = max(
                                            spath[i][j],
                                            min(spath[i][k], spath[k][j])
                                         )
    return spath


def optimal_ranking(df, d_results, borda=False):
    """
    Generate the optimized ranking by the weights calculated by the optimizer
    :param optimized_pickles: path of the outputs of the optimizer
    :param d_results: dictionary of results of the individual methods
    :param borda: whether to include the borda ranking and score in the output
    :param output_report: location of the output report
    :return: None
    """

    dict_optimal_weights = {}
    for method in d_results:

        if method in df.columns.values:
            dict_optimal_weights[method] = df[method].values[0]

    print(dict_optimal_weights)
    ranking1 = combination_ranking(d_results, dict_optimal_weights)
    df = output_to_dataframe(ranking1, d_results)

    if borda:
        d, d_scores = get_ranking_borda(d_results)
        df["RANKING_BORDA"] = df.apply(lambda row: apply_ranking_borda(d, row), axis=1)
        df["SCORE_BORDA"] = df.apply(lambda row: d_scores[row["SYMBOL"]] if row["SYMBOL"] in d_scores else np.nan,
                                     axis=1)
    return df, ranking1


def apply_ranking_borda(d, row):
    """
    :param d: dictionary of rankings
    :param row: the row of the dataframe
    :return: the ranking of the symbol
    """

    if row["SYMBOL"] in d:
        return d[row["SYMBOL"]]
    else:
        return np.nan


def get_ranking_borda(d_results):
    """
    :param d_results: Dictionary of rankings  for each individual method
    :return: dictionary of combined rankings using Borda
    """
    d_scores = {}
    d_len = {}

    # Create a dictionary of number of total elements per method
    for method in d_results.keys():
        d_len[method] = 40  # Number of genes fetch to create the pool of candidate genes

    # Now for each gene, sum the score for that gene
    for method in d_results.keys():
        for gene in d_results[method].keys():
            if not(gene in d_scores):
                d_scores[gene] = 0
            d_scores[gene] += d_len[method] - (d_results[method][gene] + 1)

    # Sort dictionary by values, the higher the better
    s_data = sorted(d_scores.items(), key=lambda item: item[1],reverse=True)
    rank, count, previous, result, score = 0, 0, None, {}, {}
    for key, num in s_data:
        count += 1
        if num != previous:
            rank += count
        previous = num
        count = 0
        result[key] = rank

    return result, d_scores
