import os
import gzip
import click
import pickle
import numpy as np
import pandas as pd
from functools import partial
from scipy.stats import dirichlet
from scipy.optimize import minimize
from scipy.optimize import basinhopping
from qc.drivers import CGC_GENES_PER_TUMOR, NEGATIVE_GENES_PER_TUMOR
from qc.deviations import Deviation
from qc.parser import Parser
from schulze import Election
from evaluation.Evaluation_Enrichment import Evaluation_Enrichment

gniter = None
gepsilon = None
gavaliable_methods = None
gdiscarded_methods = None
order_methods_ranking = ["oncodriveclust_r", "dndscv_r","oncodrivefml_r", "hotmapssignature_r","edriver_r","cbase_r","oncodriveclustl_r"]#,"mutsigcv_r"
order_methods_threshold = ["oncodriveclust_t", "dndscv_t","oncodrivefml_t", "hotmapssignature_t","edriver_t","cbase_t","oncodriveclustl_t"]#mutsigcv_t
order_methods = []
gweights = []
gt_combination = None
method_optimization = None

PATH_REGIONS = os.path.join(os.environ['SCHULZE_DATA'], '02_cds.regions.gz')
METHODS = [
    'oncodrivefml',
    'oncodriveclust',
    'dndscv',
    'hotmapssignature',
    'edriver',
    'cbase',
'oncodriveclust'

]#mutsigcv


def get_tissue(name):
    for t in NEGATIVE_GENES_PER_TUMOR.keys():
        if "_{}_".format(t) in name or name.endswith("_{}".format(t)):
            return t
    return None


class Filter:
    def __init__(self, input_run, tumor):
        self.methods = METHODS
        self.data = self.create_table(input_run, tumor)

    @staticmethod
    def statistic_outputs(input_file, tumor_name, method):
        """Calculate some statistic on the input_file
        :param input_file: path, path of the input_file
        :param tumor_name: str, name of the tumor
        :param method: str, name of the method used in the analysis
        :return: dict or None if the input file is empty
        """
        parser = Parser(method=method, gene_coordinates=PATH_REGIONS)
        df = parser.read(input_file)

        if df is None or len(df) == 0:
            return None

        positive = df[df['QVALUE'] < 0.1]
        positive = set(positive['SYMBOL']) if len(positive) > 0 else set()
        true_positive = positive & CGC_GENES_PER_TUMOR['PANCANCER']
        try:
            false_positive = positive & NEGATIVE_GENES_PER_TUMOR[get_tissue(tumor_name)]
        except KeyError as err:
            false_positive = None

        qc = Deviation(df=df, description=parser.name)
        deviation_obs_half = qc.deviation_from_null(half=True)
        areas = qc.calculate_areas()

        result = {
            'KTP': len(true_positive),
            'KFP': np.nan if false_positive is None else len(false_positive),
            'QQPLOT_STD_HALF': deviation_obs_half['deviation'],
            'QQPLOT_SLOPE_HALF': deviation_obs_half['slope'],
            'AREA_RANKING_ABSOLUTE': areas['absolute'],
            'AREA_RANKING_RELATIVE': areas['relative']
        }
        return result

    def create_table(self, input_run, tumor):
        """Run the statistics for each method on a cancer type
        :param input_run: path, location of the results folder
        :param tumor: str, name of the tumor
        :return: pandas dataframe
        """
        results = {}
        for method in self.methods:
            if method!="edriver":
                input_file = os.path.join(input_run, method, '{}.out.gz'.format(tumor))
            else:
                input_file = os.path.join(input_run, method, '{}.genes.out.gz'.format(tumor))
            res = self.statistic_outputs(input_file, tumor, method)
            if res is not None:
                results[method] = res
        return pd.DataFrame(results).T

    @staticmethod
    def find_outliers(data):
        """Find the borders of the outliers
        :param data: list of values
        :return: borders
        """
        q75, q25 = np.percentile(data, [75, 25])
        iqr = q75 - q25
        outliers_down = max(0, q25 - (iqr * 1.5))
        outliers_up = q75 + (iqr * 1.5)
        return outliers_down, outliers_up

    def filter_by(self, by):
        """Identifies methods that show abnormal values for a given metric
        :param by: str, name of the column containing the desired metric
        :return: set, discarded methods
        """
        values = self.data[by].tolist()
        methods = self.data.index.tolist()
        finite_values = [x for x in values if np.isfinite(x)]
        down, up = self.find_outliers(finite_values)
        discarded = [methods[i] for i, v in enumerate(values) if (np.isfinite(v) and v < down)]
        return set(discarded)

    def filters(self):
        """Identifies methods to be discarded by applying two filters,
        AREA_RANKING_ABSOLUTE and AREA_RANKING_RELATIVE
        :return: set, discarted methods
        """
        return self.filter_by(by='AREA_RANKING_ABSOLUTE') & self.filter_by(by='AREA_RANKING_RELATIVE')


def create_constrains():
    '''
    Create the cosntrains for the optimization process. If all methods are present, simply include the default constrains. If not, it does not include the contrain of minimum value per method.
    :param avaliable_methods:
    :return: the constrains
    '''
    # constraints of the optimization problem
    # methods=["oncodriveclust_r", "dndscv_r","oncodrivefml_r", "hotmapssignature_r","edriver_r","cbase","oncodriveclustl]
    # constraint 1: sum weights = 1;
    # constraint 2: for each weight, weight >= 0.05;
    # constraint 3: clust + hotmaps+edriver <= 0.5;
    # constraint 4: dndscv <= 0.5
    # constraint 5: fml <= 0.5

    if method_optimization == "SLSQP":
            cons = ( {'type': 'eq', 'fun': lambda w: sum(w) - 1},
                    {'type': 'ineq', 'fun': lambda w: - w[0] - w[3] -w[4] + 0.5},
                    {'type': 'ineq', 'fun': lambda w: - w[1] -w[5] + 0.5},
                    {'type': 'ineq', 'fun': lambda w: - w[2] + 0.5})
    else:
            cons = ( {'type': 'ineq', 'fun': lambda w: sum(w) - 1},
                     {'type': 'ineq', 'fun': lambda w: -sum(w) + 1},
                    {'type': 'ineq', 'fun': lambda w: - w[0] - w[3]-w[4] + 0.5},
                    {'type': 'ineq', 'fun': lambda w: - w[1] - w[5] + 0.5},
                    {'type': 'ineq', 'fun': lambda w: - w[2] + 0.5})
    if len(gavaliable_methods) == len(order_methods):

        cons = cons + ({'type': 'ineq', 'fun': lambda w: w[0] - 0.05},
        {'type': 'ineq', 'fun': lambda w: w[1] - 0.05},
        {'type': 'ineq', 'fun': lambda w: w[2] - 0.05},
        {'type': 'ineq', 'fun': lambda w: w[3] - 0.05},
        {'type': 'ineq', 'fun': lambda w: w[4] - 0.05},
        {'type': 'ineq', 'fun': lambda w: w[5] - 0.05}
        )
        return cons
    else: # Some method is missing
        list_presents = []
        for method  in order_methods:
            if method in gavaliable_methods:
                list_presents.append(method)

        l = []
        if "oncodriveclust_r" in list_presents:
            l.append( {'type': 'ineq', 'fun': lambda w: w[0] - 0.05})
        else:
            if method_optimization == "SLSQP":
                l.append({'type': 'eq', 'fun': lambda w: w[0] })
            else:
                l.append({'type': 'ineq', 'fun': lambda w: -w[0] })
                l.append({'type': 'ineq', 'fun': lambda w: w[0] })
        if "dndscv_r" in list_presents:
            l.append( {'type': 'ineq', 'fun': lambda w: w[1] - 0.05})
        else:
            if method_optimization == "SLSQP":
                l.append({'type': 'eq', 'fun': lambda w: w[1] })
            else:
                l.append({'type': 'ineq', 'fun': lambda w: -w[1] })
                l.append({'type': 'ineq', 'fun': lambda w: w[1] })
        if "oncodrivefml_r" in list_presents:
            l.append( {'type': 'ineq', 'fun': lambda w: w[2] - 0.05})
        else:
            if method_optimization == "SLSQP":
                l.append({'type': 'eq', 'fun': lambda w: w[2] })
            else:
                l.append({'type': 'ineq', 'fun': lambda w: -w[2] })
                l.append({'type': 'ineq', 'fun': lambda w: w[2] })
        if "hotmapssignature_r" in list_presents:
            l.append( {'type': 'ineq', 'fun': lambda w: w[3] - 0.05})
        else:
            if method_optimization == "SLSQP":
                l.append({'type': 'eq', 'fun': lambda w: w[3] })
            else:
                l.append({'type': 'ineq', 'fun': lambda w: -w[3] })
                l.append({'type': 'ineq', 'fun': lambda w: w[3] })

        if "edriver_r" in list_presents:
            l.append( {'type': 'ineq', 'fun': lambda w: w[4] - 0.05})
        else:
            if method_optimization == "SLSQP":
                l.append({'type': 'eq', 'fun': lambda w: w[4] })
            else:
                l.append({'type': 'ineq', 'fun': lambda w: -w[4] })
                l.append({'type': 'ineq', 'fun': lambda w: w[4] })
        if "cbase_r" in list_presents:
            l.append( {'type': 'ineq', 'fun': lambda w: w[5] - 0.05})
        else:
            if method_optimization == "SLSQP":
                l.append({'type': 'eq', 'fun': lambda w: w[5] })
            else:
                l.append({'type': 'ineq', 'fun': lambda w: -w[5] })
                l.append({'type': 'ineq', 'fun': lambda w: w[5] })

        cons = cons + tuple(l)


        return cons


def get_position(method):
    i = 0
    for method_q in order_methods:
        if method == method_q:
            return i
        i+=1


def optimize_with_seed(arg):
    '''
    Args:
        w0: array: array of weights
        func: function admitting w as argument
    Returns:
        array of weights that minimizes function
    '''

    g = arg[0]
    w = arg[1]
    cons = create_constrains()


    # hill climbing part: take w as initial value

    niter = gniter
    epsilon = gepsilon
    if method_optimization == "SLSQP":
        res = minimize(g, w, method='SLSQP', constraints=cons, options={'maxiter': niter, 'eps': epsilon,"ftol":1e-6,"disp":True})
        #minimizer_kwargs = {"constraints":cons, "method": "SLSQP"}
        #x0 = [1/6 for _ in range(6)]
        #res = basinhopping(g,x0,minimizer_kwargs=minimizer_kwargs, niter=5, T=1,stepsize=epsilon)
        #print (res)
    if method_optimization == "COBYLA":
        res = minimize(g, w, method='COBYLA', constraints=cons, options={'maxiter': niter, "tol":1e-6,"disp":True,"rhobeg":epsilon})
    return res


def optimizer(func, a=3, seeds=1, methods=['oncodrivefml', 'dndscv',  'oncodriveclust', 'oncodrivemut',"edriver","cbase","oncodriveclustl"]):
        # define symmetric dirichlet distribution with concentration parameter = a
        n = len(methods)
        alpha = [a] * n
        d = dirichlet(alpha)

        # generate an array of random seeds
        l_weights = d.rvs(size=seeds)
        # parallel hill-climbing with different seeds
        w = list(zip([func]*seeds, l_weights))
        results = []

        # FIXME Use pathos?
        # n_cores = seeds
        # with pathos.pools.ProcessPool(n_cores) as pool:
        for a in map(optimize_with_seed, w):
            results.append([a.x, a.fun])

        # sort and show the results
        print(results)
        res = sorted(results, key=lambda x: x[1], reverse=False)
        return res


def set_default_weight(d_result):
    '''
    get the function of the given metods {"oncodriveclust": w[0], "dnds": w[1], "oncodrivefml": w[2], "hotmapssignature": w[3]}
    :param d_result: the dictionary of the methods
    :return: dictionary of the default weights
    '''
    w = {}
    i = 0
    for method in d_result.keys():
        w[method] = "w["+str(i)+"]"
        i = i +1
    return w


def set_type_weights(methods,weights):
    d = {}
    i = 0
    for method in methods:
        d[method] = weights[i]
        i = i+1
    return d


def calculate_objective_function(d_ranking, weights, objective_method="Combination_Ranking",objetive_function=None,log=False):
    '''
    Given a distribution of weights returns the value of the objetive function (default: enrichment, combination_ranking)

    :param weights: dictionary of the weights where the keys are the  methods and the values are the current distribution of weights.
    :param d_ranking: dictionary of ranking of the methods where the keys are the cancer types and the methods  per cancer type.
    :param cancer: cancer type.
    :param objetive_method: objetive METHOD to optimize, default Combination_Ranking.
    :param objetive_function: function to optimize

    :return float of the current value of the objetive function for the input weights.
    '''
    global gweights
    election = Election(d_ranking)
    election.add_weights(weights)
    election.prepare()
    election.strongest_paths()
    ranking = election.combination_ranking()

    d_f = {}
    d_f[objective_method] = dict(ranking)
    if gt_combination == "RANKING":
        type = "absolute"
    else:
        type = "relative"
    d_area = objetive_function.calculate_area(d_f,type_method=type)

    return d_area[objective_method]


def prepare_output(methods_order,solutions):
    '''

    :param methods_order: methods to be output
    :param solutions: dictionary of solutions
    :return: list of solutions ready to be created by DataFrame
    '''
    l_solutions = []

    for weights, area in solutions:
        d = {}
        i = 0
        for method in methods_order:
            d[method] = weights[i]
            i = i +1
        l_solutions.append((d, area))
    return l_solutions


@click.command()
@click.option('--seeds',default=1,  help='Number of seeds.')
@click.option('--niter', default=100,help='Number of iterations')
@click.option('--epsilon', default=0.1,help='Epsilon for the optimization function')
@click.option('--foutput',type=click.Path(),help="File of output",required=True)
@click.option('--input_rankings',type=click.Path(exists=True),help="Dictionary with the ranking of the methods",required=True)
@click.option('--t_combination',help="Type of combination, ranking or threshold. Default ranking",required=True,default="RANKING")
@click.option('--optimization_algorithm',help="Algorithm for optimization. [SLSQP,COBYLA]",required=False,default="SLSQP")
@click.option('--percentage_cgc', default=1.0,help='Percentage of CGC used in the optimization. Default all.')
@click.option('--moutput',type=click.Path(exists=True),help="Input data directory", required=True)
@click.option('--cancer', type=str, required=True)
@click.option('--seed', default="T",help="Whether a seed for random generation should be defined. Default T=True. [T=True,F=False]", required=True)
def run_optimizer(seeds,niter,epsilon,foutput,input_rankings,t_combination,optimization_algorithm,percentage_cgc, moutput, cancer,seed):

    global gniter, gepsilon, gavaliable_methods, order_methods, gt_combination, method_optimization
    gniter = niter
    gepsilon = epsilon
    if seed =="T":
        np.random.seed(1)
    # Select order methods
    gt_combination = t_combination
    method_optimization = optimization_algorithm
    if t_combination == "RANKING":
        order_methods = order_methods_ranking
    else:
        order_methods = order_methods_threshold

    print(input_rankings)
    with gzip.open(input_rankings, "rb") as fd:
        d_results_methodsr = pickle.load(fd)

    print(d_results_methodsr.keys())

    # Prepare the list of available methods and the dictionary of weights
    l = []
    for method in d_results_methodsr.keys():
        l.append(method)
    gavaliable_methods = list(l)

    # Remove methods that do not reach the quality metrics
    discarded = set(["{}_r".format(m) for m in Filter(input_run=moutput, tumor=cancer).filters()])
    if len(discarded) > 0:
        print("[QC] {} discarted {}".format(cancer, discarded))

    gavaliable_methods = [m for m in gavaliable_methods if m not in discarded]

    print("Running on " + str(gavaliable_methods))
    # Set to empty disct those methods discarded or not present
    for method in order_methods:
        if not method in gavaliable_methods:
            d_results_methodsr[method] = {}

    objetive_function = Evaluation_Enrichment(percentage_cgc)

    f = partial(calculate_objective_function, d_results_methodsr, objetive_function=objetive_function)
    if t_combination == "RANKING":
        g = lambda w: -f({"oncodriveclust_r": w[0], "dndscv_r": w[1],"oncodrivefml_r": w[2], "hotmapssignature_r": w[3],"edriver_r":w[4],"cbase_r":w[5],"oncodriveclustl_r":w[6]})
    else:
        g = lambda w: -f({"oncodriveclust_t": w[0], "dndscv_r": w[1],"oncodrivefml_t": w[2], "hotmapssignature_t": w[3],"edriver_t":w[4],"cbase_t":w[5],"oncodriveclustl_t":w[6]})

    niter = gniter
    res = optimizer(g, methods=d_results_methodsr.keys(), seeds=seeds)

    output = prepare_output(order_methods, res)
    if len(output)>0:
        methods = list(output[0][0].keys())
        matrix = []
        for seed in output:
            row = []
            for method in methods:
                row.append(seed[0][method])
            row.append(seed[1])
            matrix.append(row)
        o = pd.DataFrame(matrix,columns=methods+["Objective_Function"])
    else:
        o = pd.DataFrame([],columns=order_methods+["Objective_Function"])

    o.to_csv(foutput, sep="\t", index=False, compression="gzip")


if __name__ == '__main__':

    '''
    # Example 1:
    func = lambda x: x[0] + x[2]**2
    y = func([0.2, 0.2, 0.2, 0.2, 0.2])
    x = optimize_with_seed([0.2, 0.2, 0.2, 0.2, 0.2], func)
    print(x, y)
    '''

    '''
    # Example 2:
    methods = ["oncodriveclust", "oncodriveomega", "oncodrivefml", "hotmapssignature", "mutsigcv"]
    cancer = 'LIHC'
    f = partial(calculate_objective_function, d_final, cancer_drivers, cancer)
    g = lambda w: -f({"oncodriveclust": w[0], "oncodriveomega": w[1],
                      "oncodrivefml": w[2], "hotmapssignature": w[3],
                      "mutsigcv": w[4]})
    '''
    '''
    Example 3. Run of a cancer type.

    '''

    run_optimizer()


