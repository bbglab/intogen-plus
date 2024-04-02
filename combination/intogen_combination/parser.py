
import os
from collections import defaultdict

import pandas as pd

from intogen_combination.config import CONF, METHODS


def set_ranking_genes(df_query, q_column):
    """
    Include a column with the ranking of each gene depending of the q_value.
    It allows that different genes with the same q-value are ranked in the same position
    :param df_query: Dataframe with the input
    :param q_column: name of the q-value column

    :return:
    """

    position = 0
    q_value = -1.0
    l_rankings = []
    for index, row in df_query.iterrows():
        if row[q_column] > q_value:
            position = position + 1
            q_value = row[q_column]
        l_rankings.append(position)
    df_query["Ranking"] = l_rankings
    return df_query


def create_dict_rankings(genes, rankings):
    d_out = {}
    for i in range(len(genes)):
        d_out[genes[i]] = rankings[i]
    return d_out


def parse(number_top=40, strict=True, **files):
    d = {}
    pvalues = defaultdict(dict)

    for method, file in files.items():
        c_gene, c_pvalue, c_qvalue = CONF[method]['GENE_ID'], CONF[method]['PVALUE'], CONF[method]['QVALUE']
        if os.path.exists(file):
            df = pd.read_csv(file, sep="\t")

            if df.shape[0] > 0:
                
                #FIXME No gene-name fix
                if df[c_gene].isna().any():
                    df[c_gene] = df[c_gene].fillna(df['ENSID'])
                
                assert not df[c_gene].isna().any()

                for i, r in df.iterrows():
                    try:
                        gene = r[c_gene]
                        if gene in pvalues and method in pvalues[gene]:
                            # its a repeated entry, select the one with the lowest pvalue
                            if pvalues[gene][method][0] > r[c_pvalue]:
                                pvalues[gene][method] = (r[c_pvalue], r[c_qvalue])
                        else:
                            pvalues[gene][method] = (r[c_pvalue], r[c_qvalue])
                    except KeyError as e:
                        raise e

                df = df[[c_gene, c_pvalue, c_qvalue]].drop_duplicates()
                df.sort_values(c_qvalue, inplace=True)
                df = set_ranking_genes(df, c_qvalue)
                if strict:  # do not allow draws
                    df.sort_values(c_qvalue, inplace=True)
                    df = df[(df[c_qvalue] < 1.0)].head(number_top).copy()
                else:  # include the top40 genes allowing draws
                    df = df[(df["Ranking"] < number_top) & (df[c_qvalue] < 1.0)].copy()

                genes = df[c_gene].values
                rankings = df["Ranking"].values
                d[method] = create_dict_rankings(genes, rankings)

    return d, pvalues
