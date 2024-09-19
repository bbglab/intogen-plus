import os
from collections import defaultdict
from typing import Dict, Any

import pandas as pd

from intogen_combination.config import CONF


def set_ranking_genes(df_query: pd.DataFrame, q_column: str) -> pd.DataFrame:
    """
    Include a column with the ranking of each gene depending on the q_value.
    Allows different genes with the same q-value to be ranked in the same position.

    :param df_query: DataFrame with the input
    :param q_column: Name of the q-value column
    :return: DataFrame with added 'Ranking' column
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


def create_dict_rankings(genes: list, rankings: list) -> Dict[str, int]:
    """
    Create a dictionary with gene names as keys and their rankings as values.

    :param genes: List of gene names
    :param rankings: List of corresponding rankings
    :return: Dictionary with gene rankings
    """
    d_out = {}
    for i in range(len(genes)):
        d_out[genes[i]] = rankings[i]
    return d_out


def parse(number_top: int = 40, strict: bool = True, **files: Any) -> tuple:
    """
    Parse input files to extract top-ranked genes based on q-values.

    :param number_top: Number of top-ranked genes to select
    :param strict: If True, do not allow draws (i.e., genes with the same q-value)
    :param files: Dictionary of method names and corresponding file paths
    :return: Tuple containing a dictionary of top-ranked genes for each method and a dictionary of p-values
    """
    d = {}
    pvalues = defaultdict(dict)

    for method, file in files.items():
        c_gene, c_pvalue, c_qvalue, c_ensid = (
            CONF[method]["GENE_ID"],
            CONF[method]["PVALUE"],
            CONF[method]["QVALUE"],
            CONF[method].get("ENSEMBL_ID", None)
        )

        if os.path.exists(file):
            df = pd.read_csv(file, sep="\t")

            if df.shape[0] > 0:
                # Use Ensembl ID if no gene-name where applicable
                if df[c_gene].isna().any() and c_ensid:
                        df[c_gene] = df[c_gene].fillna(df[c_ensid])

                assert not df[c_gene].isna().any()

                for i, r in df.iterrows():
                    try:
                        gene = r[c_gene]
                        if gene in pvalues and method in pvalues[gene]:
                            # It's a repeated entry, select the one with the lowest p-value
                            if pvalues[gene][method][0] > r[c_pvalue]:
                                pvalues[gene][method] = (r[c_pvalue], r[c_qvalue])
                        else:
                            pvalues[gene][method] = (r[c_pvalue], r[c_qvalue])
                    except KeyError as e:
                        raise e

                df = df[[c_gene, c_pvalue, c_qvalue]].drop_duplicates()
                df.sort_values(c_qvalue, inplace=True)
                df = set_ranking_genes(df, c_qvalue)
                if strict:  # Do not allow draws
                    df.sort_values(c_qvalue, inplace=True)
                    df = df[(df[c_qvalue] < 1.0)].head(number_top).copy()
                else:  # Include the top 'number_top' genes allowing draws
                    df = df[(df["Ranking"] < number_top) & (df[c_qvalue] < 1.0)].copy()

                genes = df[c_gene].values
                rankings = df["Ranking"].values
                d[method] = create_dict_rankings(genes, rankings)

    return d, pvalues
