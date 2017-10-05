import pandas as pd
import numpy as np
import schulze
import pickle

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
    if len(methods)==0:
        print (gene)
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

def output_to_dataframe(ranking, d_results,cancer):
    '''
    Args:
        ranking: ranking dict: dict mapping candidates to ranks
        d_results: dict mapping methods to a ranking dict
    Returns:
        dataframe encoding summary information
    '''

    l_info = []
    cgc = pd.read_csv("/workspace/datasets/CGC/generated_data2/CGCMay17_cancer_types_TCGA.tsv", sep="\t")
    cancer_drivers = cgc['Gene Symbol'].unique()
    d = d_results[cancer]
    for gene, rk in ranking.items():
        methods, best_methods, best_rank, ranks = get_voters(gene, d)
        try:
            median_rank = np.median(ranks)
            best_rank = min(ranks)
        except:
            median_rank = None
            best_rank = None
        l_info.append(
            [cancer, gene, rk, (gene in cancer_drivers), median_rank, best_rank, len(methods),
             ",".join(best_methods), ",".join(methods)])

    df_info = pd.DataFrame(l_info, columns=["Cancer_Type", "SYMBOL", "RANKING", "CGC", "Median_Ranking", "Best_Ranking",
                                            "Total_Bidders", "Highest_Bidder", "All_Bidders"])
    return df_info

if __name__ == '__main__':

    '''
    # Example 1:
    d_results = pickle.load(open("/workspace/projects/intogen/intogen4/scripts/data/dict_parsed_methods_ranking.pickle", "rb"))
    print('d_results', d_results['ESCA'].keys())
    df = output_to_dataframe(d_out, d_results)
    '''

    # Example 2:
    d_results = pickle.load(open("/workspace/projects/intogen/intogen4/scripts/data/dict_parsed_methods_threshold.pickle", "rb"))

    l_data = []
    for cancer in d_results.keys():
        election = schulze.Election(d_results[cancer])
        # election.add_weights(weights)
        election.prepare()
        election.strongest_paths_multithread(n_cores=2)
        ranking = election.combination_ranking()
        print('d_results:', d_results[cancer])
        print('ranking:', ranking)
        df = output_to_dataframe(ranking, d_results,cancer)
        df.sort_values("RANKING").to_csv("/workspace/projects/intogen/intogen4/scripts/data/results/reports/combination/threshold/"+cancer+".tsv",sep="\t",index=False)
        l_data.append(df)
    df = pd.concat(l_data)
    df.sort_values("RANKING").to_csv("/workspace/projects/intogen/intogen4/scripts/data/results/reports/combination/threshold/"+"TOTAL"+".tsv",sep="\t",index=False)

