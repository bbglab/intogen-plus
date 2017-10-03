"""This module performs simulated mutations on protein structures following the signatures associated with that cancer type and
then evaluates the clustering to detect hotspot residues.
"""
import numpy as np
import scipy.stats as stats
import src.mutations as muts
import src.pdb_structure as pstruct

NUM_MODEL_DIFF = 0
NUM_CHAIN_DIFF = 0
STRUCT_MODEL_DIFF = []
STRUCT_CHAIN_DIFF = []

def compute_pvals(density, sim_null):
    """Given the observed density, get the corresponding p-value from the
    simulated distribution.

    Parameters
    ----------
    density : list of lists, [[resNum, densityVal]]
        contains the observed density at each residue
    sim_null : np.array
        Number of instances (second column) of a specific density (first col)

    Returns
    -------
    pvals : list
        list of p-values for observed densities
    """

    # construct the cumulative distribution for the simulations
    sim_null_pval = np.cumsum(sim_null[:,1][::-1] / float(np.sum(sim_null[:,1])))[::-1]

    # get indices of p-value in simulated null distribution
    ixs = []
    for res, x in density:
        tmp_ix = 0
        for i, count in enumerate(sim_null[:,0]):
            if count <= x:
                tmp_ix = i
        ixs.append(tmp_ix)

    # number of simulations
    num_sim = float(np.sum(sim_null[:,1]))

    # look up the p-values
    #pvals = [sim_null_pval[ix] if ix < sim_null_pval.size else 1./num_sim
    pvals = [sim_null_pval[ix] if ix < sim_null_pval.size else 0.0
             for ix in ixs]

    # cdf
    sim_null[:,1] = sim_null_pval[::-1]

    return pvals, sim_null


def generate_null_dist(struct_id, model_info, chain_info,
                       cog,
                       num_mutations,
                       num_sims,
                       seed,
                       neighbours,
                       stop_criterion=np.inf,
                       max_obs=np.inf):
    """Generate a null distribution assuming a uniform rate of mutation across
    the gene.

    Parameters
    ----------
    model_info : list
        list of all possible models in the structure
    chain_info : dict
        keys are chain descriptions and values are lists of letters
        corresponding to that chain description
    cog
        center of geometry
    num_mutations : int
        number of mutations within a gene
    num_sims : int
        number of times to repeat an entire gene simulation
    seed : int
        pseudo-random number generator seed
    neighbours: dict
        dictionary of residue ids and list of
        neigbours' residue ids

    Returns
    -------
    sim_null_dist : np.array
        frequency of particular clustering observed in mutations
    """
    global NUM_MODEL_DIFF
    global NUM_CHAIN_DIFF
    global STRUCT_MODEL_DIFF
    global STRUCT_CHAIN_DIFF

    # list of all cluster numbers encountered in the
    # simulation
    sim_density = []
    prng = np.random.RandomState(seed=seed)
    res_keys = cog.keys()
    model_diff = False
    chain_diff = False

    # produce num_mutations number of mutations
    # num_sims number of times
    num_sim_gteq = 0
    for i in range(num_sims):
        tmp_mut_counts = {k: 0 for k in res_keys}

        # select a position to mutate at random
        mutated_pos_vec = prng.choice(range(len(res_keys)), size=num_mutations)
        for mutated_pos in mutated_pos_vec:
            # find equivalent chain letters to the one
            # we are mutating
            original_model = res_keys[mutated_pos][1]
            original_letter = res_keys[mutated_pos][2]
            chain_letters = pstruct.find_eq_letters(chain_info, res_keys[mutated_pos][2])
            # mutate all models with equivalent chain letters
            for m in model_info:
                for l in chain_letters:
                    position = res_keys[mutated_pos]
                    position = list(position)
                    position[1] = m
                    position[2] = l
                    position = tuple(position)
                    # check if this particular residue exists
                    # on these equivalent chains before adding
                    # a mutation
                    if position in res_keys:
                        tmp_mut_counts[position] += 1
                    elif not position[1] == original_model:
                        model_diff = True
                    elif not position[2] == original_letter:
                        chain_diff = True
                        #print("Equivalent residue does not exist on alternate chain")

        # get mutation density and add to overall sim densities
        # only for positions that have a mutation
        obs_mut_counts = {k: tmp_mut_counts[k]
                          for k in tmp_mut_counts
                          if tmp_mut_counts[k] > 0}
        mut_density = muts.mutation_density(obs_mut_counts, neighbours)
        tmp_density = [mut_density[k] for k in obs_mut_counts]
        sim_density.extend(tmp_density)

        # decide if stop early
        num_sim_gteq += sum(t>=max_obs for t in tmp_density)
        if num_sim_gteq >= stop_criterion:
            break

    if model_diff:
        NUM_MODEL_DIFF += 1
        STRUCT_MODEL_DIFF.append(struct_id)
    if chain_diff:
        NUM_CHAIN_DIFF += 1
        STRUCT_CHAIN_DIFF.append(struct_id)

    # get the frequency of each clustering pattern
    sim_null_dist = stats.itemfreq(sim_density)

    return sim_null_dist


def compute_significant_count(sim_null, sig_level):
    """
    Computes the mutation count required for clustering to be significant for
    a given radius.

    Parameters
    ----------
    sim_null : np.array
        Number of instances (second column) of a specific desnsity (first col)
    sig_level : float
        Significance level alpha required for a mutation count to be significant

    Returns
    -------
    mutation_count : int
        Number of mutations for clustering to be significant
    """

    # construct the cumulative distribution for the simulations
    sim_null_pval = np.cumsum(sim_null[:,1][::-1] / float(np.sum(sim_null[:,1])))[::-1]

    # find counts that are significant
    significant_counts = [sim_null[i][0]
                          for i in range(sim_null_pval.shape[0])
                          if sim_null_pval[i] <= sig_level]

    # return first one found if there is any
    if not significant_counts:
        # check if it was empty
        if not sim_null.any():
            return "No Mutations", 1
        # if none were found return the max count + 1
        return max([sim_null[i][0] for i in range(sim_null_pval.shape[0])]) + 1, 0
    return significant_counts[0], 1
