import Utils
import Parser
from evaluation.Evaluation_Enrichment import Evaluation_Enrichment
from bgplots.bio.enrichment import create_table
from bgplots.plot import save, show
from bgplots.matplotlib.plots import enrichment as mpl_enrichment
from bgplots.bokeh.plots import enrichment as bkh_enrichment

from tqdm import tqdm
import matplotlib.pyplot as plt
import pandas as pd
import dill as pickle
from os.path import join
import json

df_cgc = pd.read_csv("/workspace/datasets/CGC/generated_data2/CGCMay17_cancer_types_TCGA.tsv", sep="\t")
cgc = set(df_cgc['Gene Symbol'].unique())

def prep_enrichment_tables(d_results, tumor_type):
    '''
    Args:
        d_results: dict mapping cancer types into dict mapping methods into ranking dicts.
    Returns:
        d_tables: dict mapping cancer types into enrichment tables
    '''
    d_table = {}
    rankings = {}
    e = Evaluation_Enrichment(percentage=100)
    for method in d_results[tumor_type]:
        rankings[method] = e.get_list_ranking(d_results[tumor_type][method])
    d_table[tumor_type] = create_table(rankings, set_of_genes=cgc, step=1, max=100)
    return d_table

def plot_enrichment(d_results, tumor_type, label=None, plotting='mpl', folder=None):
    '''
    Args:
        d_results: dict mapping cancer types into dict mapping methods into ranking dicts
        cancer_types: iterable: all the cancer types for which plots will be produced
        plotting: str: plotting option
        path: if path is None, it displays all the plots; if path is not None, then it saves the plots in folder
    Returns:
        creates and saves one enrichment plot per cancer type
    '''
    d_tables = prep_enrichment_tables(d_results, tumor_type)
    if folder is not None:
        if label is not None:
            path_to_plot = join(folder, '{0}_{1}.png'.format(tumor_type, label))
        else:
            path_to_plot = join(folder, '{0}.png'.format(tumor_type))
    else:
        path_to_plot = None
    # title
    if label is not None:
        title = '{0}: enrichment {1}'.format(tumor_type, label)
    else:
        title = '{0}: enrichment'.format(tumor_type)
    # plotting type
    try:
        if plotting == 'mpl':
            fig = mpl_enrichment(d_tables[tumor_type], title=title)
            save(fig, path_to_plot, dpi=200, bbox_inches='tight')
        elif plotting == ' bkh':
            bkh_enrichment(d_tables[tumor_type], title=title)
    except:
        pass

if __name__ == '__main__':
    '''
    # Example 1: basic functions

    cgc = Evaluation_Enrichment.load_cgc()
    d_results_methodsr = pickle.load(
        open("/workspace/projects/intogen/intogen4/scripts/data/dict_parsed_methods_ranking.pickle", "rb"))
    d_results_methodst = pickle.load(
        open("/workspace/projects/intogen/intogen4/scripts/data/dict_parsed_methods_threshold.pickle", "rb"))
    d_results_cranking = pickle.load(
        open("/workspace/projects/intogen/intogen4/scripts/data/ranking_combination_default.pickle", "rb"))
    d_results = Utils.join_dictionaries([d_results_methodsr, d_results_methodst, d_results_cranking])
    # cancer_types = ['LIHC']
    cancer_types = None
    folder = '/workspace/projects/intogen/intogen4/scripts/data/results/plots'
    plot_enrichment(d_results, cancer_types=cancer_types, folder=folder)
    '''

    '''
    # Example 2: generate all the enrichment plots

    cgc = Evaluation_Enrichment.load_cgc()
    d_res_methodsr = pickle.load(
        open("/workspace/projects/intogen/intogen4/scripts/data/dict_parsed_methods_ranking.pickle", "rb"))
    d_res_methodst = pickle.load(
        open("/workspace/projects/intogen/intogen4/scripts/data/dict_parsed_methods_threshold.pickle", "rb"))
    d_res_cranking = pickle.load(
        open("/workspace/projects/intogen/intogen4/scripts/data/ranking_combination_default.pickle", "rb"))
    d_res_optimized = pickle.load(
        open("/workspace/projects/intogen/intogen4/scripts/data/ranking_combination_optimized.pickle", "rb"))
    d_res_cthresh = pickle.load(
        open("/workspace/projects/intogen/intogen4/scripts/data/threshold_combination_default.pickle", "rb"))
    d_res_optimized_thresh = pickle.load(
        open("/workspace/projects/intogen/intogen4/scripts/data/threshold_combination_optimized.pickle", "rb"))
    d_res = Utils.join_dictionaries([d_res_methodsr, d_res_methodst, d_res_cranking, d_res_optimized, d_res_optimized_thresh])
    d_res_rank = Utils.join_dictionaries([d_res_cranking, d_res_optimized])
    d_res_thresh = Utils.join_dictionaries([d_res_cthresh, d_res_optimized_thresh])
    cancer_types = None
    folder = '/workspace/projects/intogen/intogen4/scripts/data/results/plots'
    plot_enrichment(d_res, cancer_types=cancer_types, folder=folder)
    plot_enrichment(d_res_rank, label='rank', cancer_types=cancer_types, folder=folder)
    plot_enrichment(d_res_thresh, label='thresh', cancer_types=cancer_types, folder=folder)
    '''

    # Example 3: Compute All Enrichment Plots from Ranking Dict

    INTOGEN_PATH = '/workspace/projects/intogen/intogen4/'
    PATH_ENRICHMENT_PLOTS = join(INTOGEN_PATH, 'scripts', 'data', 'results', 'plots')
    # default set of tumor_types
    tumor_types = Parser.Parser.cancer_types
    # load the ranking dict from json
    path_to_dict = join(PATH_ENRICHMENT_PLOTS, 'combinations', 'ranking_dict.json')
    with open(path_to_dict, 'rt') as f:
        ranking_dict = json.load(f)
    # load the cgc gene set
    cgc = Evaluation_Enrichment.load_cgc()
    # loop through tumor_types
    for tumor_type in tqdm(tumor_types):
        try:
            folder = join(PATH_ENRICHMENT_PLOTS, 'combinations')
            plot_enrichment(ranking_dict, tumor_type, plotting='mpl', folder=folder)
        except:
            pass