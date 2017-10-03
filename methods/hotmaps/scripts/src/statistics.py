import numpy as np
from scipy.stats import binom

# import modules needed for logging
import logging
import os

logger = logging.getLogger(__name__)  # module logger

def cummin(x):
    """A python implementation of the cummin function in R"""
    for i in range(1, len(x)):
        if x[i-1] < x[i]:
            x[i] = x[i-1]
    return x


def bh_fdr(pval):
    """A python implementation of the Benjamani-Hochberg FDR method.

    This code should always give precisely the same answer as using
    p.adjust(pval, method="BH") in R.

    Parameters
    ----------
    pval : list or array
        list/array of p-values

    Returns
    -------
    pval_adj : np.array
        adjusted p-values according the benjamani-hochberg method
    """
    pval_array = np.array(pval)
    sorted_order = np.argsort(pval_array)
    original_order = np.argsort(sorted_order)
    pval_array = pval_array[sorted_order]

    # calculate the needed alpha
    n = len(pval)
    print n
    pval_adj = np.zeros(n)
    i = np.arange(1, n+1, dtype=float)[::-1]  # largest to smallest
    pval_adj = np.minimum(1, cummin(n/i * pval_array[::-1]))[::-1]
    return pval_adj[original_order]


def frequency_test(mut_of_interest,
                   total_mut,
                   residues_of_interest,
                   residues_at_risk):
    """Perform a binomial test on the frequency of missense mutations within
    given pre-defined residues within the gene.

    Parameters
    ----------
    mut_of_interest : {list, np.array}
        number of mutations that are deemed "of interest"
    total_mut : {list, np.array}
        total number of mutations
    residues_of_interest : {list, np.array}
        contains the number of residues of interest for a mutation.
    residues_at_risk : {list, np.array}
        contains the number of residues at risk for a mutation.

    Returns
    -------
    p_values : np.array
        p-value for each gene for binomial test
    """
    # initialize input
    p_values = np.zeros(len(mut_of_interest))
    mut = np.asarray(mut_of_interest)
    N = np.asarray(total_mut)
    residues_of_interest = np.asarray(residues_of_interest)
    residues_at_risk = np.asarray(residues_at_risk, dtype=float)
    residues_at_risk[residues_at_risk==0] = np.nan  # fill zeros to avoid divide by zero

    # calculate the background probability of mutation occurring at
    # the residues of interest
    P = residues_of_interest.astype(float) / residues_at_risk

    # iterate through each gene to calculate p-value
    logger.info('Calculating binomial test p-values . . .')
    for k in range(len(mut)):
        if not np.isnan(P[k]):
            p_val = binomial_test(mut[k], N[k], P[k])
        else:
            # catch case for nan element
            p_val = 1.0
        p_values[k] = p_val
    logger.info('Finished calculating binomial test p-values.')

    return p_values


def binomial_test(n, N, P):
    """Perform binomial test on the observed n being higher than expected.
    Specifically, N residues are at risk and of those there are n mutations
    occurred at the Np residues of interest. Given the background probability of
    a mutation at a specific residue, the p-value is calculated as the probability
    of observing n or greater mutations. Since N is large and n is small,
    it is computationally more efficient to take 1 - Pr(i<=n-1).

    Parameters
    ----------
    n : int
        number of observed mutations
    N : int
        number of residues at risk
    P : float
        background probability that a mutation would occur at a single residue

    Returns
    -------
    pval : np.array
        p-value for binomial test
    """
    if n <= 0:
        return 1.0
    pval = binom.sf(n-1, N, P)
    return pval
