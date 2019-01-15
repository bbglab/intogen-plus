import gzip
import pickle

import summary
import pandas as pd
import click

from schulze_election import combination_ranking


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


def read_optimized_dicts(optimized_pickles, d_results, output_file, output_dict, borda=False):
    """
    Generate the optmized ranking by the weights calculated by the opmitzer
    :param optimized_pickles: path of the outputs of the optimizer
    :param d_results: dictionary of results of the individual methods
    :param output_dir: directory of the output of the reports of individual cohorts
    :param output_dict: location of the output directory
    :param borda: whether to include the borda ranking and score in the output
    :param output_report: location of the output report
    :return: None
    """

    df = pd.read_csv(optimized_pickles, sep="\t", compression="gzip")

    dict_optimal_weights = {}
    for method in d_results:

        if method in df.columns.values:

            dict_optimal_weights[method] =  df[method].values[0]

    print(dict_optimal_weights)
    ranking1 = combination_ranking(d_results, dict_optimal_weights)
    df = summary.output_to_dataframe(ranking1, d_results)
    if borda:
        d, d_scores = get_ranking_borda(d_results)
        df["RANKING_BORDA"] = df.apply(lambda row: apply_ranking_borda(d, row), axis=1)
        df["SCORE_BORDA"] = df.apply(lambda row: d_scores[row["SYMBOL"]] if row["SYMBOL"] in d_scores else np.nan,
                                     axis=1)
    df.sort_values("RANKING").to_csv(output_file, sep="\t", index=False, compression="gzip")

    with gzip.open(output_dict, "wb") as fd:
        pickle.dump(ranking1, fd)

def apply_ranking_borda(d,row):
    '''

    :param d: dictionary of rankings
    :param row: the row of the dataframe
    :return: the ranking of the symbol
    '''

    if row["SYMBOL"] in d:
        return d[row["SYMBOL"]]
    else:
        return np.nan

def get_ranking_borda(d_results):
    '''

    :param d_results: Dictionary of rankings  for each individual method
    :return: dictionary of combined rankings using Borda
    '''
    d_scores = {}
    d_len = {}
    # Create a dictionary of number of total elements per method
    for method in d_results.keys():
        #d_len[method]= len(d_results[method].keys()) # This is not fair, it penalizes methods with lower number of candidate genes
        d_len[method] = 40 # Number of genes fetch to create the pool of candidate genes
    # Now for each gene, sum the score for that gene
    for method in d_results.keys():
        for gene in d_results[method].keys():
            if not(gene in d_scores):
                d_scores[gene] = 0
            d_scores[gene] += d_len[method] - (d_results[method][gene]+1)

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



def run_default_weights(d_results, output_dir, output_dict):
    '''
    Generate the optmized ranking by the weights calculated by the opmitzer
    :param dir_optimized_pickles: path of the outputs of the optimizer
    :param d_results: dictionary of results of the individual methods
    :param output_dir: directory of the output of the reports of individual cohorts
    :param output_dict: location of the output directory
    :param output_report: location of the output report
    :return: None
    '''

    num_methods = float(len(d_results.keys()))
    dict_optimal_weights = {}

    for method in d_results:
        dict_optimal_weights[method] =  1.0 / num_methods

    print(dict_optimal_weights)
    ranking1 = combination_ranking(d_results, dict_optimal_weights)

    df = summary.output_to_dataframe(ranking1, d_results)
    df.sort_values("RANKING").to_csv(output_dir, sep="\t", index=False, compression="gzip")

    with gzip.open( output_dict, "wb") as fd:
        pickle.dump(ranking1, fd)


def read_optimized_dicts_cv(dir_optimized_pickles, d_results, output_dir, output_dict):
    """
    Generate the optmized ranking by the weights calculated by the opmitzer
    :param dir_optimized_pickles: path of the outputs of the optimizer
    :param d_results: dictionary of results of the individual methods
    :param output_dir: directory of the output of the reports of individual cohorts
    :param output_dict: location of the output directory
    :param output_report: location of the output report
    :return: None
    """

    df = pd.read_csv(dir_optimized_pickles, sep="\t")
    df.sort_values("Objective_Function",inplace=True)
    dict_optimal_weights = {}
    for method in d_results:
        if method in df.columns.values:
            dict_optimal_weights[method] = df[method].values[0]

    print(dict_optimal_weights)
    ranking1 = combination_ranking(d_results, dict_optimal_weights)
    df = summary.output_to_dataframe(ranking1, d_results)
    df.sort_values("RANKING").to_csv(output_dir, sep="\t", index=False, compression="gzip")

    with gzip.open(output_dict, "wb") as fd:
        pickle.dump(ranking1, fd)


@click.command()
@click.option('--input_data',type=click.Path(exists=True),help="Dictionary with the ranking of the individual methods", required=True)
@click.option('--optimize_weights', type=click.Path(), help="Optimize weights file")
@click.option('--report_output', type=click.Path(), help="Output reports file",required=True)
@click.option('--dict_output', type=click.Path(),help="Output dictionary file", required=True)
@click.option('--type_run',help="Type of run. Optimization of the weights or default wegihts. [default,optimization,cross_validation]", required=True, default="optimization")
@click.option('--borda',help="Whether to include borda ranking in the output. Default False.", required=False, default=False)
def run_schulze(input_data, optimize_weights, report_output, dict_output, type_run, borda):

    with gzip.open(input_data, "rb") as fd:
        d_results = pickle.load(fd)

    if type_run == "optimization":
        read_optimized_dicts(optimize_weights, d_results, report_output, dict_output, borda)

    elif type_run == "cross_validation":
        read_optimized_dicts_cv(optimize_weights, d_results, report_output, dict_output)

    elif type_run == "default":
        run_default_weights(d_results, report_output, dict_output)
        return  # TODO


if __name__ == '__main__':
    run_schulze()
