from scipy.stats import dirichlet, beta
from scipy.optimize import minimize, basinhopping
import pathos.pools
from functools import partial
import pickle
from schulze import Election
from evaluation.Evaluation_Enrichment import Evaluation_Enrichment
import click
import pandas as pd
gniter = None
gepsilon = None
gavaliable_methods = None
gdiscarded_methods = None
order_methods_ranking = ["oncodriveclust_r", "oncodriveomega_r","oncodrivefml_r", "hotmapssignature_r","mutsigcv_r"]
order_methods_threshold = ["oncodriveclust_t", "oncodriveomega_t","oncodrivefml_t", "hotmapssignature_t","mutsigcv_t"]
order_methods = []
gweights = []
gt_combination = None
method_optimization = None
def create_constrains():
    '''
    Create the cosntrains for the optimization process. If all methods are present, simply include the default constrains. If not, it does not include the contrain of minimum value per method.
    :param avaliable_methods:
    :return: the constrains
    '''
    # constraints of the optimization problem
    # methods=["oncodriveclust", "oncodriveomega", "oncodrivefml", "hotmapssignature", "mutsigcv"]
    # constraint 1: sum weights = 1;
    # constraint 2: for each weight, weight >= 0.05;
    # constraint 3: mutsigcv + omega <= 0.5;
    # constraint 4: clust + hotmaps <= 0.5;
    # constraint 5: fml <= 0.5
    if method_optimization == "SLSQP":
            cons = ( {'type': 'eq', 'fun': lambda w: sum(w) - 1},
                    {'type': 'ineq', 'fun': lambda w: - w[0] - w[3] + 0.5},
                    {'type': 'ineq', 'fun': lambda w: - w[1] - w[4] + 0.5},
                    {'type': 'ineq', 'fun': lambda w: - w[2] + 0.5})
    else:
            cons = ( {'type': 'ineq', 'fun': lambda w: sum(w) - 1},
                     {'type': 'ineq', 'fun': lambda w: -sum(w) + 1},
                    {'type': 'ineq', 'fun': lambda w: - w[0] - w[3] + 0.5},
                    {'type': 'ineq', 'fun': lambda w: - w[1] - w[4] + 0.5},
                    {'type': 'ineq', 'fun': lambda w: - w[2] + 0.5})
    if len(gavaliable_methods) == len(order_methods):

        cons = cons + ({'type': 'ineq', 'fun': lambda w: w[0] - 0.05},
        {'type': 'ineq', 'fun': lambda w: w[1] - 0.05},
        {'type': 'ineq', 'fun': lambda w: w[2] - 0.05},
        {'type': 'ineq', 'fun': lambda w: w[3] - 0.05},
        {'type': 'ineq', 'fun': lambda w: w[4] - 0.05})
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
        if "oncodriveomega_r" in list_presents:
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
        if "mutsigcv_r" in list_presents:
            l.append( {'type': 'ineq', 'fun': lambda w: w[4] - 0.05})
        else:
            if method_optimization == "SLSQP":
                l.append({'type': 'eq', 'fun': lambda w: w[4] })
            else:
                l.append({'type': 'ineq', 'fun': lambda w: -w[4] })
                l.append({'type': 'ineq', 'fun': lambda w: w[4] })

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
    if method_optimization == "COBYLA":
        res = minimize(g, w, method='COBYLA', constraints=cons, options={'maxiter': niter, "tol":1e-6,"disp":True,"rhobeg":epsilon})
    return res

def optimizer(func, a=3, seeds=1,methods=['oncodrivefml', 'oncodriveomega', 'mutsigcv', 'oncodriveclust', 'oncodrivemut']):
        # define symmetric dirichlet distribution with concentration parameter = a
        n = len(methods)
        alpha = [a] * n
        d = dirichlet(alpha)

        # generate an array of random seeds
        l_weights = d.rvs(size=seeds)
        # parallel hill-climbing with different seeds
        w = list(zip([func]*seeds, l_weights))
        results = []
        n_cores = seeds
        with pathos.pools.ProcessPool(n_cores) as pool:
            for a in pool.imap(optimize_with_seed, w):
                results.append([a.x, a.fun])

        # sort and show the results
        print (results)
        res = sorted(results, key=lambda x: x[1], reverse=False)
        return res


def set_default_weight(d_result):
    '''
    get the function of the given metods {"oncodriveclust": w[0], "oncodriveomega": w[1], "oncodrivefml": w[2], "hotmapssignature": w[3],"mutsigcv": w[4]}
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

def calculate_objective_function(d_ranking, cancer, weights, objective_method="Combination_Ranking",objetive_function=None,  methods=["oncodriveclust","oncodriveomega","oncodrivefml","hotmapssignature","mutsigcv"],log=False):
    '''
    Given a distribution of weights returns the value of the objetive function (default: enrichment, combination_ranking)

    :param weights: dictionary of the weights where the keys are the  methods and the values are the current distribution of weights.
    :param d_ranking: dictionary of ranking of the methods where the keys are the cancer types and the methods  per cancer type.
    :param cancer: cancer type.
    :param objetive_method: objetive METHOD to optimize, default Combination_Ranking.
    :param objetive_function: function to optimize
    :param  methods:  list of methods to optimize the weights.
    :return float of the current value of the objetive function for the input weights.
    '''
    global gweights
    d_query = {}
    d_query[cancer] = {}
    election = Election(d_ranking[cancer])
    #weightsn = set_type_weights(methods,weights)

    election.add_weights(weights)
    election.prepare()
    #election.strongest_paths_multithread(n_cores=1)
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

    for weights,area in solutions:
        d = {}
        i = 0
        for method in methods_order:
            d[method] = weights[i]
            i = i +1
        l_solutions.append((d,area))
    return l_solutions

def read_discarded_methods(file_discarded,ttype,t_combination):
    '''
    Read the discarded methods by quality control and remove then from the voting system
    :param file_discarded: File with the tabulated file with the information
    :param ttype: tumor type
    :param t_combination: type of combination. RANKING or THRESHOLD
    :return: the name of discarded methods, if any
    '''

    df = pd.read_csv(file_discarded,sep="\t")

    methods = df[(df["TUMOR"]==ttype)&(df["METRIC"].str.contains(t_combination))]["METHOD"].values
    l = list()
    for method in methods:
        if t_combination == "RANKING":
            l.append(method+"_r")
        else:
            l.append(method+"_t")
    return l


@click.command()
@click.option('--seeds',default=1,  help='Number of seeds.')
@click.option('--niter', default=100,help='Number of iterations')
@click.option('--epsilon', default=0.1,help='Epsilon for the optimization function')
@click.option('--foutput',type=click.Path(),help="File of output",required=True)
@click.option('--cancer',help="Cancer type to run the input",required=True)
@click.option('--input_rankings',type=click.Path(exists=True),help="Dictionary with the ranking of the methods",required=True)
@click.option('--discarded_methods',type=click.Path(),help="File with the discarded methods for each cancer type",required=False)
@click.option('--t_combination',help="Type of combination, ranking or threshold. Default ranking",required=True,default="RANKING")
@click.option('--optimization_algorithm',help="Algorithm for optimization. [SLSQP,COBYLA]",required=False,default="SLSQP")
@click.option('--percentage_cgc', default=1.0,help='Percentage of CGC used in the optimization. Default all.')
#@click.option('--log',type=click.Path(),help="Log file with the weights of the optimization steps",required=False,default=None)
#/workspace/projects/intogen/intogen4/runs/intogen4_20170614/plots/discarted_analyses.txt
def run_optimizer(seeds,niter,epsilon,foutput,cancer,input_rankings,discarded_methods,t_combination,optimization_algorithm,percentage_cgc):

    global gniter,gepsilon,gavaliable_methods,order_methods,gt_combination,method_optimization
    gniter = niter
    gepsilon = epsilon
    print ("Tumor type",cancer)
    print ("Directory output",foutput)
    # Select order methods
    gt_combination = t_combination
    method_optimization = optimization_algorithm
    if t_combination == "RANKING":
        order_methods = order_methods_ranking
    else:
        order_methods = order_methods_threshold
    d_results_methodsr= pickle.load( open( input_rankings, "rb" ) )

    # Prepare the list of available methods and the dictionary of weights
    l = []
    for method in d_results_methodsr[cancer].keys():
        l.append(method)
    gavaliable_methods = list(l)
    # Remove methods that do not reach the quality metrics
    discarded = read_discarded_methods(discarded_methods,cancer,t_combination)
    for method in gavaliable_methods:
        if method in discarded:
            gavaliable_methods.remove(method)
    print ("Running on " + str(gavaliable_methods))
    # Set to empty disct those methods discarded or not present
    for method in order_methods:
        if not method in gavaliable_methods:
            d_results_methodsr[cancer][method] = {}

    objetive_function = Evaluation_Enrichment(percentage_cgc)

    f = partial(calculate_objective_function, d_results_methodsr, cancer,objetive_function=objetive_function,methods=list(d_results_methodsr[cancer].keys()))
    if t_combination == "RANKING":
        g = lambda w: -f({"oncodriveclust_r": w[0], "oncodriveomega_r": w[1],"oncodrivefml_r": w[2], "hotmapssignature_r": w[3],"mutsigcv_r": w[4]})
    else:
        g = lambda w: -f({"oncodriveclust_t": w[0], "oncodriveomega_t": w[1],"oncodrivefml_t": w[2], "hotmapssignature_t": w[3],"mutsigcv_t": w[4]})

    niter = gniter
    res = optimizer(g,methods=d_results_methodsr[cancer].keys(),seeds = seeds)
    output = prepare_output(order_methods,res)
    if len(output)>0:
        methods = list(output[0][0].keys())
        matrix = []
        for seed in output:
            row = []
            for method in methods:
                row.append(seed[0][method])
            row.append(seed[1])
            matrix.append(row)
        o =  pd.DataFrame(matrix,columns=methods+["Objective_Function"])
    else:
        o =  pd.DataFrame([],columns=order_methods+["Objective_Function"])
    #pickle.dump(output,open(foutput,"wb"))

    #o.to_pickle(foutput+cancer+"_"+t_combination+".pickle")
    o.to_csv(foutput+cancer+"_weights"+".tsv",sep="\t",index=False)


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


