
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import uniform, chi2, combine_pvalues
from statsmodels.stats.multitest import multipletests

from intogen_combination.config import METHODS


def trunc(pvals, epsilon=1e-10):
    """trunc to prevent arbitrarily low pvals"""

    mask = (pvals < epsilon)
    pvals[mask] = epsilon
    return pvals


def impute(pvals):
    """impute with uniform [0,1] distribution"""

    mask1 = np.isnan(pvals)
    mask2 = (pvals == 1)
    mask = mask1 | mask2
    pvals[mask] = uniform.rvs(size=len(pvals[mask]))
    return pvals


def distance_matrix(df, metric='correlation'):
    """
    df has columns of the form 'PVALUE_<method>'
    returns the matrix of distance between methods
    """

    methods_list = METHODS
    x = df[['PVALUE_' + m for m in methods_list]].values
    x = impute(x)
    Y = pdist(x.T, metric=metric)
    return squareform(Y, force='no', checks=True)


def fisher(pvals, var=None):
    return combine_pvalues(pvals, method='fisher')[1]


def brown(pvalues, var=None):
    """
    :param pvalues: array of pvalues
    :param var: parameter for Brown's method
    :return: combined p-value
    """

    k = len(pvalues)
    e = 2 * k
    f = (e ** 2) / var
    c = k / f
    psi = -2 * sum(np.log(pvalues))
    return 1 - chi2.cdf(psi / c, 2 * f)


def custom_combination(df, comb):
    """
    :param df: summary table with results of all the methods
    :param comb: combination method: 'fisher' or 'brown'
    :return summary table with new combination columns
    Remark: KIRC data returned good clustering of driver discovery methods
    """

    pval_cols = list(map(lambda x: '_'.join(['PVALUE', x]), METHODS))
    D = distance_matrix(df[pval_cols])

    # compute parameters to feed Brown's method
    k = len(pval_cols)
    cov = np.zeros((k, k))
    for i, c in enumerate(pval_cols):
        for j in range(5):
            if i < j:
                c = 1 - D[i, j]
                cov[i, j] = (3.263 * c) + (0.710 * c ** 2) + (0.027 * c ** 3)
    var = (4 * k) + (2 * np.sum(cov))

    g = globals()  # required to get functions by name from the global namespace

    df['PVALUE_' + comb] = df[['PVALUE_' + m for m in METHODS]].apply(lambda x: g[comb](impute(x), var=var), axis=1)
    df['PVALUE_trunc_' + comb] = df[['PVALUE_' + m for m in METHODS]].apply(lambda x: g[comb](trunc(impute(x)), var=var), axis=1)
    df['QVALUE_' + comb] = multipletests(df['PVALUE_' + comb].values, method='fdr_bh')[1]
    df['QVALUE_trunc_' + comb] = multipletests(df['PVALUE_trunc_' + comb].values, method='fdr_bh')[1]
    return df
