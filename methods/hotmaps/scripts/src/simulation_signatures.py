"""This module performs simulated mutations on mutations and
then evaluates the clustering to detect hotspot residues.
"""
import numpy as np
import scipy.stats as stats
import src.mutations as muts
import src.pdb_structure as pstruct
import randomizer_aa
from src import simulate_mutations_signatures as sig
import pandas as pd
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


# TODO Profile this function
def generate_null_dist_sig(struct_id,coordinates,model_info, chain_info,
                       cog,
                       num_mutations,
                       num_sims,
                       seed,
                       neighbours,signatures,cancer_type,d_correspondence,
                       stop_criterion=np.inf,
                       max_obs=np.inf):
    """Generate a null distribution assuming a uniform rate of mutation across
    the gene.

    Parameters
    ----------
    coordinates: list of dictionaries of the coordinates of the query protein
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
    signature: 
    path to the signature file
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
   
    res_keys = cog.keys()
    model_diff = False
    chain_diff = False

    # produce num_mutations number of mutations
    # num_sims number of times
    num_sim_gteq = 0
    
    
    # select a position to mutate at random

    results = randomizer_aa.randomize_region(number_mutations=num_mutations, input_regions=coordinates, number_simulations=num_sims,signatures=signatures,cancer_type=cancer_type,cores=4)
    #print results[0]
    mutated_pos_vec = sig.map_generated_mutations(struct_id, results,d_correspondence)
    if len(mutated_pos_vec) == 0:
        return []
    
    for key in mutated_pos_vec.keys():
        tmp_mut_counts = {k: 0 for k in res_keys}
        for (chain,mutated_pos) in mutated_pos_vec[key]: # tru of 
            # find equivalent chain letters to the one

            original_model = 0
            original_letter = chain
            chain_letters = pstruct.find_eq_letters(chain_info, chain)
            # mutate all models with equivalent chain letters
            for m in model_info:
                for l in chain_letters:
                    position = (struct_id, m, l, (' ',mutated_pos,' '))
                    # check if this particular residue exists
                    # on these equivalent chains before adding
                    # a mutation
                    if position in res_keys:
                        tmp_mut_counts[position] += 1
                    elif not position[1] == original_model:
                        model_diff = True
                    elif not position[2] == original_letter:
                        chain_diff = True

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
        
    if model_diff:
        NUM_MODEL_DIFF += 1
        STRUCT_MODEL_DIFF.append(struct_id)
        
    if chain_diff:
        NUM_CHAIN_DIFF += 1
        STRUCT_CHAIN_DIFF.append(struct_id)

    # get the frequency of each clustering pattern
    sim_null_dist = stats.itemfreq(sim_density)

    return sim_null_dist


def read_file_coordinates(file_cordinates):
    """
    Read the file of genomic coordinates of each protein ID and returns a DataFrame with the information

    Parameters
    ----------
    file_coordinates: path to the file of coordinates
    Returns
    -------
    DataFrame with the genomic coordinates indexed by ProteinName
    """
    df = pd.read_csv(file_cordinates,names=["chromosome","start","end","strand","pdb_id","chain"],sep="\t")
    #df = pd.read_csv(file_cordinates,names=["chromosome","start","end","strand","geneid","category"],sep="\t")
    #df.set_index("geneid",inplace=True)
    return df
    

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
