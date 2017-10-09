import gzip
from collections import defaultdict
from functools import partial
import operator
#import pathos.pools
import pickle
import glob
import re
import summary
import pandas as pd
import click

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

class Election(object):
    """
    dict mapping voter to dict mapping candidates to valid ranks
    ties are allowed
    """
    def __init__(self, ballot_dict):
        self.dict = ballot_dict
        self.all_candidates = list(set.union(*[set(self.dict[voter].keys()) for voter in self.dict]))
        self.pref = defaultdict(lambda: defaultdict(float))
        self.spath = defaultdict(lambda: defaultdict(float))
        self.weights = None

    def get_voter_list(self):
        voter_list = []
        for ballot in self.dict:
            voter_list.append(ballot.get_voter())
        return voter_list

    def get_ballot(self, voter):
        return Ballot({voter: self.dict[voter]})

    def add_ballot(self, ballot):
        self.dict.update(ballot.dict)

    def add_weights(self, weight_dict):
        self.weights = weight_dict

    def prepare(self):
        """
        Returns:
            self.pref: dict mapping candidates to dict mapping candidates to weighted
            votes supporting primary key has higher rank than secondary key.
        """
        if self.weights is None:
            weights = {}
            for voter in self.dict:
                weights[voter] = 1.
            self.weights = dict(weights)

        for voter in self.dict:
            d = self.dict[voter]
            for i in self.all_candidates:
                if i not in d.keys():
                    for j in d:
                        self.pref[j][i] += self.weights[voter]
                else:
                    r = d[i]
                    for j in d:
                        if d[j] < r:
                            self.pref[j][i] += self.weights[voter]

    def strongest_paths(self):
        """
        Returns:
            dict mapping candidates to dict mapping candidates to highest strength of path
            joining primary key and secondary key.

        Notice that this requires prior run of method self.prepare()
        """
        for i in self.all_candidates:
            for j in self.all_candidates:
                if i != j:
                    if self.pref[i][j] > self.pref[j][i]:
                        self.spath[i][j] = self.pref[i][j]
        for i in self.all_candidates:
            for j in self.all_candidates:
                if i != j:
                    for k in sorted(self.all_candidates):
                        if (i != k) and (j != k):
                            self.spath[j][k] = max(
                                                    self.spath[j][k],
                                                    min(self.spath[j][i], self.spath[i][k])
                                                   )

    def strongest_paths_multithread(self, n_cores=1):
        """
        Multiprocessing version of self.strongest_paths() method -- see above.
        n_cores = number of cores used in multiprocessing
        """
        for i in self.all_candidates:
            for j in self.all_candidates:
                if i != j:
                    if self.pref[i][j] > self.pref[j][i]:
                        self.spath[i][j] = self.pref[i][j]
        f = partial(strongest_paths_by_chunk, self.all_candidates, self.spath)
        n_chunks = n_cores
        chunk_list = chunkizate(self.all_candidates, n_chunks)
        results = defaultdict(lambda: defaultdict(float))

        # FIXME Use pools ??
        #with pathos.pools.ProcessPool(n_cores) as pool:
        for a in map(f, chunk_list):
            results.update(a)
        self.spath = results

    def scores_dict(self):
        """
        For each candidate C, the "score" of C is defined as the number of paths
        from any other candidate to C that are stronger than its reverse path.

        Returns:
            dict mapping candidates to scores

        Notice that this method requires prior run of method self.strongest_paths()
        """
        scores_dict = {}
        for i in self.all_candidates:
            score = 0
            for j in self.spath[i]:
                if self.spath[i][j] < self.spath[j][i]:
                    score += 1
            scores_dict[i] = score
        return sorted(scores_dict.items(), key=operator.itemgetter(1), reverse=True)

    def combination_ranking(self):
        """
        Returns:
            dict mapping candidate with its rank after combination;
            notice that ties are allowed in the resulting ranking.
        """
        sorted_scores = self.scores_dict()
        ranking = {}
        prev_score = None
        prev_rank = None
        counter = 1
        while len(sorted_scores) > 0:
            c = sorted_scores.pop()
            if prev_score is None:
                ranking[c[0]] = 1
                prev_score = c[1]
                prev_rank = 1
            elif prev_score == c[1]:
                ranking[c[0]] = prev_rank
            elif prev_score < c[1]:
                ranking[c[0]] = counter
                prev_score = c[1]
                prev_rank = counter
            counter += 1

        return ranking

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


def read_optimized_dicts(optimized_pickles, d_results, output_file, output_dict):
    """
    Generate the optmized ranking by the weights calculated by the opmitzer
    :param optimized_pickles: path of the outputs of the optimizer
    :param d_results: dictionary of results of the individual methods
    :param output_dir: directory of the output of the reports of individual cohorts
    :param output_dict: location of the output directory
    :param output_report: location of the output report
    :return: None
    """

    df = pd.read_csv(optimized_pickles, sep="\t", compression="gzip")

    dict_optimal_weights = {}
    for method in d_results:

        if method in df.columns.values:

            dict_optimal_weights[method] =  df[method].values[0]

    election = Election(d_results)
    election.add_weights(dict_optimal_weights)

    print ( dict_optimal_weights)
    election.prepare()
    election.strongest_paths()
    ranking1 = election.combination_ranking()
    df = summary.output_to_dataframe(ranking1, d_results, cancer_type)

    df.sort_values("RANKING").to_csv(output_file, sep="\t", index=False, compression="gzip")

    with gzip.open(output_dict, "wb") as fd:
        pickle.dump(ranking1, fd)


def run_default_weights(d_results,output_dir,output_dict,name_run,suffix):
    '''
    Generate the optmized ranking by the weights calculated by the opmitzer
    :param dir_optimized_pickles: path of the outputs of the optimizer
    :param d_results: dictionary of results of the individual methods
    :param output_dir: directory of the output of the reports of individual cohorts
    :param output_dict: location of the output directory
    :param output_report: location of the output report
    :return: None
    '''
    d_total = {}
    l_data = []

    for cancer_type in d_results.keys():

            d_total[cancer_type] = {}
            num_methods = float(len(d_results[cancer_type].keys()))
            dict_optimal_weights = {}
            for method in d_results[cancer_type]:

                dict_optimal_weights[method] =  1.0 / num_methods

            election = Election(d_results[cancer_type])
            election.add_weights(dict_optimal_weights)
            print (cancer_type)
            print ( dict_optimal_weights)
            election.prepare()
            election.strongest_paths()
            ranking1 = election.combination_ranking()
            d_total[cancer_type][name_run] = ranking1#d_total[cancer_type]["Combination_Threshold_Optimized"] = ranking1
            df = summary.output_to_dataframe(ranking1, d_results,cancer_type)
            df.sort_values("RANKING").to_csv(output_dir+"/"+cancer_type+".tsv",sep="\t",index=False)#df.sort_values("RANKING").to_csv("/workspace/projects/intogen/intogen4/scripts/data/results/reports/optimization/threshold/"+cancer_type+".tsv",sep="\t",index=False)
            l_data.append(df)
    df = pd.concat(l_data)
    df.sort_values("RANKING").to_csv(output_dir+"/"+"TOTAL"+".tsv",sep="\t",index=False)#df.sort_values("RANKING").to_csv("/workspace/projects/intogen/intogen4/scripts/data/results/reports/optimization/threshold/"+"TOTAL"+".tsv",sep="\t",index=False)
    pickle.dump( d_total, open( output_dict, "wb" ) )#pickle.dump( d_total, open( "/workspace/projects/intogen/intogen4/scripts/data/threshold_combination_optimized.pickle", "wb" ) )


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

    election = Election(d_results)
    election.add_weights(dict_optimal_weights)
    print(dict_optimal_weights)
    election.prepare()
    election.strongest_paths()
    ranking1 = election.combination_ranking()
    df = summary.output_to_dataframe(ranking1, d_results, cancer_type)
    df.sort_values("RANKING").to_csv(output_dir, sep="\t", index=False, compression="gzip")

    with gzip.open(output_dict, "wb") as fd:
        pickle.dump( ranking1, fd)


@click.command()
@click.option('--input_data',type=click.Path(exists=True),help="Dictionary with the ranking of the individual methods", required=True)
@click.option('--optimize_weights', type=click.Path(), help="Optimize weights file")
@click.option('--report_output', type=click.Path(), help="Output reports file",required=True)
@click.option('--dict_output', type=click.Path(),help="Output dictionary file", required=True)
@click.option('--type_run',help="Type of run. Optimization of the weights or default wegihts. [default,optimization,cross_validation]",required=True,default="optimization") # Example "Combination_Ranking_Optimized"
@click.option('--type_input',help="Type of input. Threshold or Ranking [ranking,threshold]",required=True, default="ranking") # Example "Combination_Ranking_Optimized"
def run_schulze(input_data, optimize_weights, report_output, dict_output, type_run, type_input):

    with open(input_data, "rb") as fd:
        d_results = pickle.load(fd)

    if type_run == "optimization":
        read_optimized_dicts(optimize_weights, d_results, report_output, dict_output)

    elif type_run == "cross_validation":
        read_optimized_dicts_cv(optimize_weights, d_results, report_output, dict_output)

    elif type_run == "default":
        run_default_weights(d_results, report_output, dict_output, name, type_input)
        return  # TODO





if __name__ == '__main__':
    '''
    Example 1. Toy model.
    x = Ballot({'mutsigcv': {'A': 1, 'B': 2, 'C': 3, 'D': 4}})
    y = Ballot({'oncodrivefml': {'A': 1, 'B': 2, 'C': 2, 'D': 4}})
    z = Ballot({'omega1': {'A': 2, 'B': 2, 'C': 2, 'D': 1}})
    x.validate()
    y.validate()
    z.validate()
    election = Election(x.dict)
    election.add_ballot(y)
    election.add_ballot(z)
    weights = {'mutsigcv': 0, 'oncodrivefml': 0, 'omega1': 1}
    election.add_weights(weights)
    election.prepare()
    print(election.dict)
    # election.strongest_paths()
    election.strongest_paths_multithread(n_cores=2)
    print(election.spath)
    ranking = election.combination_ranking()
    print(ranking)
    print(chunkizate([1,2,3,4,5,6,7,8,9,10,11], 4))
    '''

    '''
    Example 2.LIHC, [mutsigcv,oncodrivefml, oncodriveomega]


    d_results= pickle.load( open( "/workspace/projects/intogen/intogen4/scripts/data/dict_parsed_methods_ranking.pickle", "rb" ) )
    weights = {'mutsigcv_r': 0.182, 'oncodrivefml_r': 0.185, 'oncodriveomega_r': 0.20,"oncodriveclust_r":0.12,"hotmapssignature_r":0.30} # Default
    weights ={'hotmapssignature_r': 0.099315862543989927, 'oncodrivefml_r': 0.13913990142371621, 'mutsigcv_r': 0.24596849605476062, 'oncodriveomega_r': 0.32088383525594577, 'oncodriveclust_r': 0.19469190472158751}
    cancers = ["LIHC"]

    for cancer in cancers:
        election = Election(d_results[cancer])
        election.add_weights(weights)
        election.prepare()
        print(election.dict)
        # election.strongest_paths_multithread(n_cores=2)
        election.strongest_paths()
        print(election.spath)
        ranking1 = election.combination_ranking()
        ranking2 = election.combination_ranking()
    #pickle.dump( ranking1, open( "/workspace/projects/intogen/intogen4/scripts/data/test_3.pickle", "wb" ) )
    print([x[0] for x in sorted(ranking1.items(),key=lambda x: (x[1],x[0]))])
    e = Evaluation_Enrichment()
    dict_auc1 = e.calculate_area({"Combination_Ranking": ranking1})

    print (dict_auc1)
    '''
    '''
    Run for all cancer types with default weights


    d_results= pickle.load( open( "/workspace/projects/intogen/intogen4/scripts/data/dict_parsed_methods_threshold_simulated.pickle", "rb" ) )


    cancer_types = ["ACC","BLCA","BRCA","CESC","CHOL","COAD","COADREAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"]
    d_ranking = {}
    for cancer in cancer_types:
        print (cancer)
        d_ranking[cancer] = {}
        election = Election(d_results[cancer])

        election.prepare()
        election.strongest_paths_multithread(n_cores=4)
        ranking = election.combination_ranking()
        d_ranking[cancer]["Combination_Threshold"] = ranking

    pickle.dump( d_ranking, open( "/workspace/projects/intogen/intogen4/scripts/data/threshold_combination_simulated.pickle", "wb" ) )
    '''
    '''
    Read the optimized weights and generate the runs
    python /workspace/projects/intogen/intogen4/scripts/Schulze/schulze.py --input_data /workspace/projects/intogen/intogen4/scripts/data/dict_parsed_methods_threshold.pickle --directory_optimize_weights /workspace/projects/intogen/intogen4/scripts/data/results/optimization/weights/ --directory_output  /workspace/projects/intogen/intogen4/scripts/data/results/optimization/threshold --dict_output /workspace/projects/intogen/intogen4/scripts/data/results/optimization/threshold/dict_threshold_combination_optimized.pickle --name Combination_Threshold_Optimized --type_run optimization --type_input threshold

    '''
    run_schulze()
