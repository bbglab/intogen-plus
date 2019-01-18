import csv
import gzip
import os
from collections import defaultdict

import pandas as pd
import pickle
import click


class Parser():

    methods = ["hotmaps","oncodrivefml","dndscv","cbase","oncodriveclustl","smregions"]
    column_keys = {
        "hotmaps":     ["GENE",        "q-value",      "Min p-value"],
        "oncodrivefml":         ["SYMBOL",      "Q_VALUE",      "P_VALUE"],
        "dndscv":               ["gene_name",   "qallsubs_cv",  "pallsubs_cv"],
        "cbase":                ["gene",        "q_pos",    "p_pos"],
        "oncodriveclustl":      ["SYMBOL",         "Q_ANALYTICAL",       "P_ANALYTICAL"],
        "smregions": ["HUGO_SYMBOL", "Q_VALUE", "P_VALUE"]

    }

    def __init__(self, path, thresholds={"hotmaps":0.1,"oncodrivefml":0.1,"dndscv":0.1,"cbase":0.1,"oncodriveclustl":0.1,"smregions":0.1} ):
        self.path = path
        self.thresholds = thresholds

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

    def create_dictionary_outputs(self, number_top=40, strict=True, cancer=None):
        '''

        :param number_top: Number of top selected ranking genes. (default: top 40).
        :param strict: Whether the number_top is strict or relative to the ranking. (default: strict)
        :param cancer: cancer type
        :return: a dictionary with the the structure {"Cancer":{"Method":{"Gene1":1,"Gene2":2}}
        '''

        d = {}
        pvalues = defaultdict(dict)
        for method in self.methods:
            path = os.path.join(self.path, method, "{}.out.gz".format(cancer))
            
            if os.path.exists(path):
                df = pd.read_csv(path, sep="\t")

                if df.shape[0]>0:

                    for i, r in df.iterrows():
                        try:
                            pvalues[r[self.column_keys[method][0]]][method] = (r[self.column_keys[method][2]],r[self.column_keys[method][1]])
                        except KeyError as e:
                            raise e

                    df = df[self.column_keys[method]].drop_duplicates()
                    q_value_c = self.column_keys[method][1] # Q-value column name
                    genes_c = self.column_keys[method][0] # gene name column name
                    df.sort_values(q_value_c,inplace=True)
                    df = self.set_ranking_genes(df,q_value_c)
                    if not strict: # include the top40 genes allowing draws
                        df = df[(df["Ranking"]<number_top)&(df[q_value_c]<1.0)].copy()
                    else: # do not allow draws
                        df.sort_values(q_value_c,inplace=True)
                        df = df[(df[q_value_c]<1.0)].head(number_top).copy()
                    genes = df[genes_c].values
                    rankings = df["Ranking"].values
                    d[method] = Parser.create_dict_rankings(genes,rankings)

        return d, pvalues

    @staticmethod
    def get_value(dict_values, key, position=0):

        if key in dict_values:
            return dict_values[key][position] # Position = 0 for pvalue, Position = 1 for Qvalue
        else:
            return None


@click.command()
@click.option('--input',type=click.Path(exists=True),help="Input data directory", required=True)
@click.option('--cancer', type=str, required=True)
@click.option('--output', type=click.Path(),help="Output file", required=True)
def run_parser(input, cancer, output):

    p = Parser(input)
    d_outr, pvalues = p.create_dictionary_outputs(cancer=cancer)
    with gzip.open("{}".format(output), "wb") as fd:
        pickle.dump(d_outr, fd)
    print (d_outr)
    with gzip.open("{}b".format(output), "wt") as fd:
        writer = csv.writer(fd, delimiter='\t')
        writer.writerow(['SYMBOL'] + ["PVALUE_{}".format(h) for h in p.methods]+["QVALUE_{}".format(h) for h in p.methods])

        for gene, values in pvalues.items():
            writer.writerow([gene] + [Parser.get_value(values,m,0) for m in p.methods] + [Parser.get_value(values,m,1) for m in p.methods] )


if __name__ == "__main__":
    run_parser()

