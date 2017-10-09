import gzip
import os
import pandas as pd
import pickle
import click


class Parser():

    methods = ["hotmapssignature","oncodrivefml","mutsigcv","oncodriveomega","oncodriveclust"] # "oncodriveclust","oncodriveomega"
    column_keys = {"hotmapssignature":["GENE","q-value"],"oncodrivefml":["SYMBOL","Q_VALUE"],"mutsigcv":["gene","q"],"oncodriveomega":["SYMBOL","q_value"],"oncodriveclust":["SYMBOL","QVALUE"]}

    def __init__(self, path, thresholds ={"hotmapssignature":0.1,"oncodrivefml":0.1,"mutsigcv":0.1,"oncodriveomega":0.1,"oncodriveclust":0.1} ):
        self.path = path
        self.thresholds = thresholds

    @staticmethod
    def read_hugo(path=os.environ['SCHULZE_DATA']):
        with open(os.path.join(path, "ENSEMBL_HUGO_1807107.pickle"),"rb") as fd:
            return pickle.load(fd)

    @staticmethod
    def create_dict_rankings(genes,rankings):
        d_out = {}
        for i in range(len(genes)):
            d_out[genes[i]] = rankings[i]
        return d_out

    def set_ranking_genes(self, df_query,q_column):
        '''
        Include a column with the ranking of each gene depending of the q_value. It allows that different genes with the same q-value are ranked in the same position
        :param df: Dataframe with the input
        :param q_column: name of the q-value column

        :return:
        '''

        position = 0
        q_value = -1.0
        l_rankings = []
        for index,row in df_query.iterrows():
            if row[q_column] > q_value:
                position = position + 1
                q_value = row[q_column]
            l_rankings.append(position)
        df_query["Ranking"] = l_rankings
        return df_query

    def create_dictionary_outputs(self, type_selection="threshold", number_top=40, strict=True, cancer=None):
        '''

        :param type_selection: Type of selection of the genes  [threshold,ranking] (default: threshold).
        :param number_top: Number of top selected ranking genes. (default: top 40).
        :param strict: Whether the number_top is strict or relative to the ranking. (default: strict)
        :return: a dictionary with the the structure {"Cancer":{"Method":{"Gene1":1,"Gene2":2}}
        '''
        d_hugo = Parser.read_hugo()
        suffix = type_selection[0]

        d = {}
        for method in self.methods:
            path = os.path.join(self.path, method, "{}.out.gz".format(cancer))
            if os.path.exists(path):
                df = pd.read_csv(path, sep="\t")
                if df.shape[0]>0:
                    if method == "oncodriveomega" or method == "oncodriveclust":
                        # Include the Hugo_Symbol
                        df["SYMBOL"] = df.apply(lambda row: str(d_hugo[row["GENE"]]) if row["GENE"] in d_hugo else "-" ,axis=1 )

                    df = df[self.column_keys[method]].drop_duplicates()
                    q_value_c = self.column_keys[method][1] # Q-value column name
                    genes_c = self.column_keys[method][0] # gene name column name
                    df.sort_values(q_value_c,inplace=True)
                    df = self.set_ranking_genes(df,q_value_c)
                    if type_selection=="threshold": # Select genes below the threshold
                        df = df[df[q_value_c]<=self.thresholds[method]].copy()
                        if method == "oncodriveomega":
                            df = df.head(limit_mutsigcv).copy()
                        if method == "oncodriveclust":
                            df = df.head(limit_hts).copy()
                        if method == "mutsigcv":
                            limit_mutsigcv = df.shape[0]
                        if method == "hotmapssignature":
                            limit_hts = df.shape[0]
                    else: # Select the ranking genes. In case there is a tie, returns the top 40 + the tied ones.
                        if strict == False:
                            df = df[(df["Ranking"]<number_top)&(df[q_value_c]<1.0)].copy()
                        else:
                            df.sort_values(q_value_c,inplace=True)
                            df = df[(df[q_value_c]<1.0)].head(number_top).copy()#strict

                    genes = df[genes_c].values
                    rankings = df["Ranking"].values
                    d[method+"_"+suffix] = Parser.create_dict_rankings(genes,rankings)

        return d


@click.command()
@click.option('--input',type=click.Path(exists=True),help="Input data directory", required=True)
@click.option('--cancer', type=str, required=True)
@click.option('--output', type=click.Path(),help="Output file", required=True)
@click.option('--selection', help="[ranking,threshold]", default="ranking")
def run_parser(input, cancer, output, selection):
    p = Parser(input)
    d_outr = p.create_dictionary_outputs(type_selection=selection, cancer=cancer)
    with gzip.open("{}".format(output), "wb") as fd:
        pickle.dump(d_outr, fd)


if __name__ == "__main__":
    run_parser()
    '''
    INTOGEN RUN

    p = Parser("/projects_bg/bg/shared/projects/intogen/intogen4/runs/intogen4_20170614/")
    d_outr = p.create_dictionary_outputs(type_selection="ranking")
    pickle.dump( d_outr, open( "/workspace/projects/intogen/intogen4/scripts/data/dict_parsed_methods_ranking.pickle", "wb" ) )
    d_outt = p.create_dictionary_outputs(type_selection="threshold")
    pickle.dump( d_outt, open( "/workspace/projects/intogen/intogen4/scripts/data/dict_parsed_methods_threshold.pickle", "wb" ) )
    '''

    '''
    SIMULATED

    p = Parser("/projects_bg/bg/shared/projects/intogen/intogen4/runs/simulated_20170621/")
    d_results_r = p.create_dictionary_outputs(type_selection="ranking")
    pickle.dump( d_results_r, open( "/workspace/projects/intogen/intogen4/scripts/data/dict_parsed_methods_ranking_simulated.pickle", "wb" ) )
    d_results_t = p.create_dictionary_outputs(type_selection="threshold")
    pickle.dump( d_results_t, open( "/workspace/projects/intogen/intogen4/scripts/data/dict_parsed_methods_threshold_simulated.pickle", "wb" ) )
    p = Parser("/projects_bg/bg/shared/projects/intogen/intogen4/runs/simulated_20170621/",thresholds ={"hotmapssignature":1.0,"oncodrivefml":1.0,"mutsigcv":1.0,"oncodriveomega":1.0,"oncodriveclust":1.0})
    d_results_total = p.create_dictionary_outputs(type_selection="threshold")
    pickle.dump( d_results_total, open( "/workspace/projects/intogen/intogen4/scripts/data/dict_parsed_methods_all_simulated.pickle", "wb" ) )
    '''


