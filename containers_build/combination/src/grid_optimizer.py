import os
import gzip
import click
import pickle
import numpy as np
import pandas as pd
import itertools
from scipy.optimize import basinhopping
from functools import partial

from qc.deviations import Deviation
from qc.parser import Parser
from schulze_election import combination_ranking
from evaluation.Evaluation_Enrichment import Evaluation_Enrichment

from configobj import ConfigObj

# Global variables
CODE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)))
configs_file = os.path.join(CODE_DIR, 'config', 'intogen_qc.cfg')
configs = ConfigObj(configs_file)
METHODS = configs["methods"].keys() # reads them in order

gavaliable_methods = None

PATH_REGIONS = os.path.join(os.environ['INTOGEN_DATASETS'], 'shared', 'cds.regions.gz')


class Filter:

    def __init__(self, input_run, tumor):
        self.methods = METHODS
        self.data = self.create_table(input_run, tumor)

    @staticmethod
    def statistic_outputs(input_file, tumor_name, method):
        """
        Calculate some statistic on the input_file
        :param input_file: path, path of the input_file
        :param tumor_name: str, name of the tumor
        :param method: str, name of the method used in the analysis
        :return: dict or None if the input file is empty
        """
        parser = Parser(method=method, gene_coordinates=PATH_REGIONS)
        df = parser.read(input_file)

        if df is None or len(df) == 0:
            return None

        qc = Deviation(df=df, description=parser.name)
        areas = qc.calculate_areas(up_to=40)

        result = {
            'AREA_RANKING_ABSOLUTE': areas['absolute'],
            'AREA_RANKING_RELATIVE': areas['relative']
        }
        return result

    def create_table(self, input_run, tumor):
        """
        Run the statistics for each method on a cancer type
        :param input_run: path, location of the results folder
        :param tumor: str, name of the tumor
        :return: pandas dataframe
        """
        results = {}
        for method in self.methods:
            input_file = os.path.join(input_run, method, '{}.out.gz'.format(tumor))            
            res = self.statistic_outputs(input_file, tumor, method)

            if res is not None:
                results[method] = res
            else:
                results[method] = {
                    'AREA_RANKING_ABSOLUTE': 0,
                    'AREA_RANKING_RELATIVE': 0
                }

        return pd.DataFrame(results).T

    @staticmethod
    def find_outliers(data):
        """
        Find the borders of the outliers
        :param data: list of values
        :return: borders
        """
        q75, q25 = tuple(map(lambda x: np.percentile(data, x), [75, 25]))
        iqr = q75 - q25
        outliers_down = max([0, q25 - (iqr * 1.5)])
        outliers_up = q75 + (iqr * 1.5)
        return outliers_down, outliers_up

    def filter_by(self, by):
        """
        Identifies methods that show abnormal values for a given metric
        :param by: str, name of the column containing the desired metric
        :return: set, discarded methods
        """
        values = self.data[by].tolist()
        methods = self.data.index.tolist()
        finite_values = [x for x in values if np.isfinite(x)]
        down, up = self.find_outliers(finite_values)
        discarded = [methods[i] for i, v in enumerate(values) if (np.isfinite(v) and v <= down)]
        return set(discarded)

    def filters(self):
        """
        Identifies methods to be discarded by applying two filters,
        AREA_RANKING_ABSOLUTE and AREA_RANKING_RELATIVE
        :return: set, discarded methods
        """
        return self.filter_by(by='AREA_RANKING_ABSOLUTE') & self.filter_by(by='AREA_RANKING_RELATIVE')


def calculate_objective_function(d_ranking, objective_function, weights, objective_method="Combination_Ranking"):
    """
    Given a distribution of weights returns the value of the objective function
    (default: enrichment, combination_ranking)

    :param weights: dictionary of the weights where the keys are the
                    methods and the values are the current distribution of weights.
    :param d_ranking: dictionary of ranking of the methods where the keys are the
                    cancer types and the methods per cancer type.
    :param objective_method: objective METHOD to optimize, default Combination_Ranking.
    :param objective_function: function to optimize

    :return float of the current value of the objective function for the input weights.
    """
    ranking = combination_ranking(d_ranking, weights)
    d_f = {objective_method: dict(ranking)}
    d_area = objective_function.calculate_area(d_f, type_method="absolute")
    return d_area[objective_method]


def prepare_output(methods_order, solutions):
    """
    :param methods_order: methods to be output
    :param solutions: dictionary of solutions
    :return: list of solutions ready to be created by DataFrame
    """
    l_solutions = []
    for weights, area in solutions:
        d = {}
        for i, method in enumerate(methods_order):
            d[method] = weights[i]
        l_solutions.append((d, area))
    return l_solutions

# constraint function wrappers


def array_component(w, i):
    return w[i]


def lower_bound(w, i):
    return w[i] - 0.05


def upper_bound(w, i):
    return - w[i] + 0.5


def simplex_bound(w):
    return sum(w) - 1


def clustering_bound(w):
    return - w[0] - w[3] - w[4] + 0.5


def recurrence_bound(w):
    return - w[1] - w[5] + 0.5


def fm_bound(w):
    return - w[2] + 0.5


def grid_optimize(func, low_quality=set()):
    """
    :param: func: function to be optimized
    :return: best candidates

    These are the constraints that must be in place:
    constraint 1: sum w_i = 1
    constraint 2: w_clust + w_hotmaps +w_smregions <= 0.5
    constraint 3: w_dndscv + w_cbase <= 0.5
    constraint 4: w_fml <= 0.5
    constraint 5*: w_i <= 0.05 for all w_i may not apply if some w_i is discarded

    """

    optimum = {k: None for k in METHODS}
    optimum.update({'Objective_Function': 0})
    low_quality_index = None
    if len(low_quality) > 0:
        low_quality_index = [METHODS.index(v) for v in low_quality]
    low = [0.05 if v not in low_quality else 0 for v in METHODS]  # change lower bound on discarded methods
    for w in itertools.product(np.linspace(0, 1, 21), repeat=5):
        if sum(w) <= 0.95:  # belongs to simplex
            w = list(np.append(w, [1 - sum(w)]))
            if low_quality_index is not None:  # add zeros at discarded / low-quality methods
                for ind in low_quality_index:
                    w[ind] = 0.
                total_weight = sum(w)
                w = [v / total_weight for v in w]
            if (low[0] <= w[0] < 0.5) and (low[3] <= w[3] < 0.5) and (low[4] <= w[4] < 0.5):
                if w[0] + w[3] + w[4] < 0.5:  # cluster constraint
                    if low[2] <= w[2] < 0.5:    # fm bias constraint
                        if (low[1] <= w[1] < 0.5) and (low[5] <= w[5] < 0.5):
                            if w[1] + w[5] < 0.5:  # recurrence constraint
                                f = func(w)
                                if optimum['Objective_Function'] > f:
                                    optimum['Objective_Function'] = f
                                    for i, v in enumerate(METHODS):
                                        optimum[v] = w[i]
    return optimum


def create_scipy_constraints(low_quality=set()):
    """
    Create the constraints for the optimization process.
    """

    # constraints of the optimization problem
    # methods=["oncodriveclustl_r", "dndscv_r", "oncodrivefml_r", "hotmaps_r", "cbase"]
    # constraint 1: sum weights = 1
    # constraint 2: for each weight, weight >= 0.05
    # constraint 3: clust + hotmaps + e_driver <= 0.5
    # constraint 4: dndscv <= 0.5
    # constraint 5: fml <= 0.5

    cons = [{'type': 'eq', 'fun': simplex_bound},
            {'type': 'ineq', 'fun': clustering_bound},
            {'type': 'ineq', 'fun': recurrence_bound},
            {'type': 'ineq', 'fun': fm_bound}]
    for ind, v in enumerate(METHODS):
        if v in low_quality:
            cons += [{'type': 'eq', 'fun': partial(array_component, i=ind)}]
        else:
            cons += [{'type': 'ineq', 'fun': partial(lower_bound, i=ind)}]
            cons += [{'type': 'ineq', 'fun': partial(upper_bound, i=ind)}]
    return tuple(cons)


def satisfy_constraints(w, low_quality=set()):
    satisfy = True
    for i, v in enumerate(METHODS):
        if v in low_quality:
            satisfy = satisfy and (abs(array_component(w, i)) < 0.01)
        else:
            satisfy = satisfy and (lower_bound(w, i) >= 0)
    satisfy = satisfy and (fm_bound(w) >= 0)
    satisfy = satisfy and (clustering_bound(w) >= 0)
    satisfy = satisfy and (recurrence_bound(w) >= 0)
    satisfy = satisfy and (abs(simplex_bound(w)) < 0.01)
    return satisfy


def optimize_with_seed(func, w0, low_quality=set()):
    """
    Args:
        func: function admitting w as argument
        w0: array: array of weights
        low_quality: list of methods for which weight shall be 0.0.
    Returns:
        array of weights that minimizes function
    """

    cons = create_scipy_constraints(low_quality=low_quality)
    niter = 25
    epsilon = 0.02
    options = {'maxiter': niter, 'eps': epsilon, 'ftol': 1e-3, 'disp': True}
    minimizer_kwargs = {'method': 'SLSQP', 'options': options, 'constraints': cons}
    res = basinhopping(func, w0, minimizer_kwargs=minimizer_kwargs, niter=1, stepsize=0.05)
    return res


def full_optimizer(cancer, input_rankings, method_reject, moutput, percentage_cgc, seed):
    '''

    :param cancer:
    :param input_rankings:
    :param method_reject:
    :param moutput:
    :param percentage_cgc:
    :param seed:
    :return:
    '''

    global gavaliable_methods
    if seed == 'T':
        np.random.seed(1)

    # Select order methods
    with gzip.open(input_rankings, "rb") as fd:
        d_results_methodsr = pickle.load(fd)
    print(d_results_methodsr.keys())

    # Prepare the list of available methods and the dictionary of weights
    l = []
    for method in d_results_methodsr.keys():
        l.append(method)
    gavaliable_methods = list(l)

    # Remove methods that do not reach the quality metrics
    if not (method_reject is None):
        discarded = set()
        discarded.add(method_reject)
    else:
        discarded = set(["{}".format(m) for m in Filter(input_run=moutput, tumor=cancer).filters()])

    # Include discarded from command line
    if len(discarded) > 0:
        print("[QC] {} discarded {}".format(cancer, discarded))
    gavaliable_methods = [m for m in gavaliable_methods if m not in discarded]
    print("Running on " + str(gavaliable_methods))

    # Set to empty those methods discarded or not present
    for method in METHODS:
        if method not in gavaliable_methods:
            d_results_methodsr[method] = {}


    # Instantiate the enrichment objective function
    objective_function = Evaluation_Enrichment(percentage_cgc)
    f = partial(calculate_objective_function, d_results_methodsr, objective_function)

    def func(w):
        return -f(dict(zip(METHODS, w)))

    # best solution in 1/20 resolution grid, augmented with basin-hopping/SLSQP optimization
    grid_optimum = grid_optimize(func, low_quality=discarded)  # get optimum candidate in the grid
    w = np.array([grid_optimum[k] for k in METHODS])
    res = optimize_with_seed(func, w)  # basin-hopping/SLSQP using grid optimum candidate as initial guess
    res_dict = dict(zip(METHODS, list(res.x)))
    res_dict['Objective_Function'] = res.fun

    # choose the best one grid or basin-hopping, unless basin-hopping does not fulfill the constraints
    if res_dict['Objective_Function'] > grid_optimum['Objective_Function']:
        out_df = pd.DataFrame({k: [v] for k, v in grid_optimum.items()})
    else:
        r = np.array([res_dict[k] for k in METHODS])
        if satisfy_constraints(r, low_quality=discarded):
            out_df = pd.DataFrame({k: [v] for k, v in res_dict.items()})
        else:
            out_df = pd.DataFrame({k: [v] for k, v in grid_optimum.items()})
    return out_df


def skip_optimizer(method_reject, moutput, cancer):

    global gavaliable_methods

    gavaliable_methods = METHODS

    # Remove methods that do not reach the quality metrics
    if not (method_reject is None):
        discarded = set()
        discarded.add(method_reject + "_r")
    else:
        discarded = set(["{}_r".format(m) for m in Filter(input_run=moutput, tumor=cancer).filters()])

    # Include discarded from command line
    gavaliable_methods = [m for m in gavaliable_methods if m not in discarded]
    print("Running on " + str(gavaliable_methods))

    # Create a uniform vector of weights
    N = len(gavaliable_methods)
    df = pd.DataFrame({k: [1/N] for k in gavaliable_methods})
    if len(discarded) > 0:
        print("[QC] {} discarded {}".format(cancer, discarded))
        for md in discarded:
            df[md] = 0.0
    df["Objective_Function"] = np.nan
    return df



@click.command()
@click.option('--foutput',
              type=click.Path(),
              help="File of output",
              required=True)
@click.option('--input_rankings',
              type=click.Path(exists=True),
              help="Dictionary with the ranking of the methods",
              required=True)
@click.option('--percentage_cgc',
              default=1.0,
              help='Percentage of CGC used in the optimization. Default all.')
@click.option('--moutput',
              type=click.Path(exists=True),
              help="Input data directory",
              required=True)
@click.option('--method_reject',
               type=click.STRING,
               help="Method for which weight in combination is forced to be zero",
               required=False)
@click.option('--cancer',
              type=str,
              required=True)
@click.option('--seed',
              default="T",
              type=str,
              help="Whether a seed for random generation should be defined. Default T=True. [T=True,F=False]",
              required=True)
def run_optimizer(foutput, input_rankings, percentage_cgc, moutput, cancer, seed, method_reject):

    if float(percentage_cgc) > 0.0:
        out_df = full_optimizer(cancer, input_rankings, method_reject, moutput, percentage_cgc, seed)
    else:
        out_df = skip_optimizer(method_reject, moutput, cancer)

    # write to table
    out_df.to_csv(foutput, sep="\t", index=False, compression="gzip")


if __name__ == '__main__':
    run_optimizer()