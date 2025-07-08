import numpy as np
from src.density import *
from src import utils
import Bio.PDB
import string
import os

def find_neighbors(cog, max_dist=10):
    """
    Find neigbours of each residue within
    a MaxDist neigbourhood

    Parameters:
    -----------
    cog : dictionary
        dictionary of residue ids and their
        centers of geomeetry
    max_dist : int
        radius of neighborhood in angstroms

    Returns:
    --------
    neigbours_dict : dictionary
        dictionary of residue ids and a list
        of their neighbors' residue ids
    """
    # For every residue, find all residues within a MaxDist
    # angstrom radius.
    res_keys = cog.keys()
    neighbors_dict = {}
    for Residue1 in res_keys:
        neighbors = []
        for Residue2 in cog.keys():
            # make sure we dont count the residue itself
            # as a neighbour
            if (Residue1 == Residue2):
                continue
            tmp_dist = distance([cog[Residue1],
                                 cog[Residue2]])
            if tmp_dist <= (max_dist):
                neighbors.append(Residue2)
        neighbors_dict[Residue1] = neighbors

    return neighbors_dict


def find_neighbors_for(cog, target_res, max_dist=10):
    """
    Overloaded version of find neighbors
    that allows us to find the neighbors
    of only the residues we want

    Parameters:
    -----------
    cog : dictionary
        dictionary of residue ids and their
        centers of geomeetry
    max_dist : int
        radius of neighborhood in angstroms
    target_res : list
        list of residues whose neighbors we're looking for
    Returns:
    --------
    neigbours_dict : dictionary
        dictionary of residue ids and a list
        of their neighbors' residue ids
    """
    # For every residue, find all residues within a MaxDist
    # angstrom radius.
    res_keys = cog.keys()
    neighbors_dict = {}
    for Residue1 in target_res:
        neighbors = []
        for Residue2 in cog.keys():
            # make sure we dont count the residue itself
            # as a neighbour
            if (Residue1 == Residue2):
                continue
            tmp_dist = distance([cog[Residue1],
                                 cog[Residue2]])
            if tmp_dist <= (max_dist):
                neighbors.append(Residue2)
        neighbors_dict[Residue1] = neighbors

    return neighbors_dict


def find_neighbors_1D(cog, max_codon_dist=3):
    """
    Find neigbours of each residue within
    a MaxDist neigbourhood in 1D space

    Parameters:
    -----------
    cog : dictionary
        dictionary of residue ids and their
        centers of geomeetry
    max_dist : int
        distance in number of residues considered a
        neighbourhood

    Returns:
    --------
    neigbours_dict : dictionary
        dictionary of residue ids and a list
        of their neighbors' residue ids
    """
    # For every residue, find all residues within a MaxDist
    # angstrom radius.
    res_keys = cog.keys()
    neighbors_dict = {}
    for Residue1 in res_keys:
        neighbors = []
        for Residue2 in cog.keys():
            # make sure we dont count the residue itself
            # as a neighbour
            if (Residue1 == Residue2):
                continue

            # figure out the codon position of the residue
            residue1_pos = int(Residue1[3][1])
            residue1_chain = Residue1[2]
            residue1_model = Residue1[1]
            residue2_pos = int(Residue2[3][1])
            residue2_chain = Residue2[2]
            residue2_model = Residue2[1]

            # how many codons are they apart
            tmp_dist = abs(residue1_pos - residue2_pos)

            if residue1_model == residue2_model and residue1_chain == residue2_chain and tmp_dist <= (max_codon_dist):
                neighbors.append(Residue2)
        neighbors_dict[Residue1] = neighbors

    return neighbors_dict


def find_neighbors_for_1D(cog, target_res, max_codon_dist=3):
    """
    Overloaded version of find neighbors 1D
    that allows us to find the neighbors
    of only the residues we want

    Parameters:
    -----------
    cog : dictionary
        dictionary of residue ids and their
        centers of geomeetry
    max_dist : int
        radius of neighborhood in angstroms
    target_res : list
        list of residues whose neighbors we're looking for
    Returns:
    --------
    neigbours_dict : dictionary
        dictionary of residue ids and a list
        of their neighbors' residue ids
    """
    # For every residue, find all residues within a MaxDist
    # angstrom radius.
    res_keys = cog.keys()
    neighbors_dict = {}
    for Residue1 in target_res:
        neighbors = []
        for Residue2 in cog.keys():
            # make sure we dont count the residue itself
            # as a neighbour
            if (Residue1 == Residue2):
                continue


            # figure out the codon position of the residue
            residue1_pos = int(Residue1[3][1])
            residue1_chain = Residue1[2]
            residue1_model = Residue1[1]
            residue2_pos = int(Residue2[3][1])
            residue2_chain = Residue2[2]
            residue2_model = Residue2[1]

            # how many codons are they apart
            tmp_dist = abs(residue1_pos - residue2_pos)

            if residue1_model == residue2_model and residue1_chain == residue2_chain and tmp_dist <= (max_codon_dist):
                neighbors.append(Residue2)

        neighbors_dict[Residue1] = neighbors

    return neighbors_dict


def fix_chain_letter(chain):
    """Fix the chain letter when no letter is assigned.

    Parameters
    ----------
    chain : Biopython chain
        chain object from biopython structure

    Returns
    -------
    chain : biopython chain
        chain with fixed letter

    """
    # if the chain id is ' '
    # we can conlcude there is only one chain
    # and its id is 'A'
    if chain.get_id() == ' ':
        chain.id = 'A'
        #chain_id = 'A'
        #full_id = list(full_id)
        #full_id[2] = 'A'
        #full_id = tuple(full_id)
    #else:
        #chain_id = chain.get_id()
    return chain


def get_structure_info(structure, mchains, mres, mcount, struct_chains, ttype_ixs):
    """
    Gets relevant structure info such as centers of geometry
    of all residues, mutated residues, and all models in the structure

    Parameters:
    ----------
    structure : BioPDB.Structure
        the structure to be parsed
    mchains : list
        chains containing mutations
    mres : list
        residues containing mutations
    mcount : list
        number of mutations present
        in mutated chains and residue pairings
    struct_chains : list
        list of chain letters in info file
        compiled from Mupit database
    ttype_ixs : dictionary
        dictionary of tumour types and indices
        that correspond to mutation information

    Returns:
    --------
    structure_info : tuple
        (mut_res_centers_of_geometry : dictionary of dictionaries
            for each tumour type the residue ids and corresponding cogs,
        mut_res_mutation_counts : dictionary of dictionaries
            for each tumour type the residue ids and correspinding mut cts,
        all_res_centers_of_geometry : dictionary,
        models : list)
    """

    # get center of geometry of each residue in the structure
    # also check if the residue has a mutation present
    # create two dictionaries with full residue id as keys
    # and mutation counts and center of geometry as values
    # note: obtain residue id with residue.id[1]
    #       or with residue.get_full_id()[3][1]
    #       in order to check whether an residue
    #       is a hetero residue, check that
    #       residue.get_full_id()[3][0] is ' '
    #       if it contains 'H_something' or 'W'
    #       we discard it as a hetero atom
    mut_res_centers_of_geometry = {k:{} for k in ttype_ixs}
    mut_res_mutation_counts = {k:{} for k in ttype_ixs}
    all_res_centers_of_geometry = {}
    models = []
    # observed_chains = []

    for model in structure:
        models.append(model.get_id())

        for chain in model:

            # this is done to prevent double counting
            # mutations that are present on different models
            # but on the same chain, and residue
            #if chain.get_id() in observed_chains:
            #    continue
            #observed_chains.append(chain.get_id())

            # sometimes homology models have invalid chain letters
            chain = fix_chain_letter(chain)
            chain_id = chain.get_id()

            for residue in chain:

                # ignore hetero residues
                if not (residue.get_full_id()[3][0] == ' '):
                    continue

                full_id = residue.get_full_id()

                # if the chain id is ' '
                # we can conlcude there is only one chain
                # and its id is 'A'
                #if chain.get_id() == ' ':
                    #chain_id = 'A'
                    #full_id = list(full_id)
                    #full_id[2] = 'A'
                    #full_id = tuple(full_id)
                #else:
                    #chain_id = chain.get_id()

                # calculate center of geometry
                center_of_geometry = np.sum(atom.coord for atom in residue) / len(residue)

                for tumour in ttype_ixs:
                    #if tumour == 'PANCAN':
                        #continue
                    for i in ttype_ixs[tumour]:
                        # check if this residue in this chain has a mutation
                        if mchains[i] == chain_id and mres[i] == full_id[3][1]:
                            num_mutations = mcount[i]
                            if full_id[2] in struct_chains:
                                if full_id in mut_res_mutation_counts[tumour]:
                                    mut_res_mutation_counts[tumour][full_id] += num_mutations
                                else:
                                    mut_res_centers_of_geometry[tumour][full_id] = center_of_geometry
                                    mut_res_mutation_counts[tumour][full_id] = num_mutations
                                    break

                # insert the information into the dictionary
                if (full_id[2] in struct_chains):
                    all_res_centers_of_geometry[full_id] = center_of_geometry

    ttypes = mut_res_mutation_counts.keys()
    return mut_res_centers_of_geometry, mut_res_mutation_counts, all_res_centers_of_geometry, models


def calc_center_of_geometry(structure, struct_chains):
    """Calculates the center of geometry."""
    all_res_centers_of_geometry = {}

    for model in structure:
        for chain in model:
            # sometimes homology models have invalid chain letters
            chain = fix_chain_letter(chain)

            for residue in chain:
                # ignore hetero residues
                if not (residue.get_full_id()[3][0] == ' '):
                    continue

                center_of_geometry = np.sum(atom.coord for atom in residue) / len(residue)

                # insert the information into the dictionary
                full_id = residue.get_full_id()
                if (full_id[2] in struct_chains):
                    all_res_centers_of_geometry[full_id] = center_of_geometry
    return all_res_centers_of_geometry


def get_filtered_atom_list(struct):
    """Filters out unwanted atoms in PDB structure.

    Removes heteroatoms from atoms when finding neighbor residues.

    Parameters
    ----------
    struct : Bio.PDB structure
        biopython structure object

    Returns
    -------
    atom_list : list
        list of atoms after filtering
    """
    atom_list = []
    for model in struct:
        for chain in model:
            for residue in chain:
                # ignore hetero residues
                if not (residue.get_full_id()[3][0] == ' '):
                    continue

                # add residue atoms
                atom_list.extend(residue.get_list())
    return atom_list


def get_interface_residues(structure, radius):
    """Find amino acid residues at the interface either between
    another protein, DNA, or RNA.

    Parameters
    ----------
    structure : Bio.PDB structure
        structure object to find interface residues for
    radius : float
        radius close enough to define contact

    Returns
    -------
    interface_res : dict
        dictionary annotating each residue as interface or not
    """
    # get resnames for nucleotides
    dna_resnames = ['DT', 'DA', 'DC', 'DG']
    rna_resnames = ['A', 'C', 'G', 'T', 'U']

    # get list of atoms without heteroatoms
    atoms = get_filtered_atom_list(structure)

    # examine neighbor pairs to find interface residues
    interface_res = {}
    neighbors = Bio.PDB.NeighborSearch(atoms)
    res_pairs = neighbors.search_all(radius, level='R')
    for res1, res2 in res_pairs:
        res1_full_id = res1.get_full_id()
        res2_full_id = res2.get_full_id()

        # skip if on same chain in same model
        if (res1_full_id[1] == res2_full_id[1]) and \
            (res1_full_id[2] == res2_full_id[2]):
            continue

        # hack to prevent multiple reporting of residues
        # on different models but the same actual chain
        #res1_full_id, res2_full_id = list(res1_full_id), list(res2_full_id)
        #res1_full_id[1], res2_full_id[1] = 0, 0
        #res1_full_id, res2_full_id = tuple(res1_full_id), tuple(res2_full_id)

        # find if interface is protein-protein, protein-dna, or protein-rna
        # convention is 1 indicates interface residue for protein-protein,
        # protein-dna, and protein-rna in the saved list, respecitively.
        is_res1_aa = Bio.PDB.is_aa(res1)
        is_res2_aa = Bio.PDB.is_aa(res2)
        if is_res1_aa and is_res2_aa:
            interface_res.setdefault(res1_full_id, [0, 0, 0])
            interface_res.setdefault(res2_full_id, [0, 0, 0])
            interface_res[res1_full_id][0] = 1
            interface_res[res2_full_id][0] = 1
        elif is_res1_aa and res2.get_resname().strip() in dna_resnames:
            interface_res.setdefault(res1_full_id, [0, 0, 0])
            interface_res[res1_full_id][1] = 1
        elif is_res2_aa and res1.get_resname().strip() in dna_resnames:
            interface_res.setdefault(res2_full_id, [0, 0, 0])
            interface_res[res2_full_id][1] = 1
        elif is_res1_aa and res2.get_resname().strip() in rna_resnames:
            interface_res.setdefault(res1_full_id, [0, 0, 0])
            interface_res[res1_full_id][2] = 1
        elif is_res2_aa and res1.get_resname().strip() in rna_resnames:
            interface_res.setdefault(res2_full_id, [0, 0, 0])
            interface_res[res2_full_id][2] = 1

    return interface_res


def get_buried_residues(structure, cutoff, tmp_dir, dssp_path):
    """Finds buried residues by using relative solvent accessible surface area.

    """
    # get structure id
    structure_id = structure.id

    all_letters = set(string.ascii_uppercase) | set(string.ascii_lowercase)

    # flatten models into a single model due to limitations of DSSP
    id_map = {}
    for k, model in enumerate(structure):
        if k == 0:
            #used_letters = set(model.child_dict.keys())
            used_letters = set()
            for chain in model:
                if chain.get_id() == ' ':
                    chain.id = 'A'
            #for l in used_letters:
                id_map[(model.id, chain.id)] = (model.id, chain.id)
                used_letters.add(chain.id)
            new_model = model.id
        else:
            for chain in model:
                left_over = all_letters - used_letters
                if not left_over:
                    # if run out of chain letters just return nothing
                    return []
                new_letter = left_over.pop()
                used_letters.add(new_letter)
                old_letter = chain.id
                chain.id = new_letter
                id_map[(new_model, new_letter)] = (model.id, old_letter)

                # add numbers if there is not more letters left
                if not (all_letters - used_letters):
                    all_letters.update(set(string.digits) | set(string.punctuation))

            model.id = new_model

    # save new structure to tmp dir
    io = Bio.PDB.PDBIO()
    io.set_structure(structure)
    tmp_path = os.path.join(tmp_dir, structure_id+'.pdb')
    io.save(tmp_path)

    # read in tmp structure
    tmp_structure = utils.read_structure(tmp_path, structure_id, quiet=True)

    # find the solvent accessibility for residues
    dssp_results = Bio.PDB.DSSP(tmp_structure[0], tmp_path, dssp=dssp_path)

    # get bfactors for each amino acid residue
    bfacs_missing = [r
                     for r in tmp_structure.get_residues()
                     if Bio.PDB.is_aa(r) and 'CA' not in r.child_dict]
    bfacs = [r['CA'].get_bfactor()
             for r in tmp_structure.get_residues()
             if Bio.PDB.is_aa(r) and 'CA' in r.child_dict]
    mean_bfac = np.mean(bfacs)
    std_bfac = np.std(bfacs)

    # format output
    output = []
    for result in dssp_results:
        # skip if not an amino acid
        if not Bio.PDB.is_aa(result[0]):
            continue

        # format the ID
        full_id = result[0].get_full_id()
        #if full_id[2] == ' ':
            #full_id[2] = 'A'
        try:
            orig_model_chain = list(id_map[full_id[1:3]])
        except:
            print(full_id, id_map)
            raise
        # fix missing letter for homology models
        if orig_model_chain[1] == ' ':
            orig_model_chain[1] = 'A'

        # record whether it was buried
        if 'CA' in result[0].child_dict:
            norm_bfactor = (result[0]['CA'].get_bfactor() - mean_bfac) / std_bfac
        else:
            norm_bfactor = None
        line = [structure_id] + orig_model_chain + [result[0].id[1], result[3], norm_bfactor]
        if result[3] <= cutoff:
            line.append(1)
        else:
            line.append(0)
        output.append(line)

    # delete tmp file
    if os.path.exists(tmp_path): os.remove(tmp_path)

    return output


def find_eq_letters(chain_info, letter):
    """
    Find chain letters with the same chain description

    Parameters
    ----------
    chain_info : dict
        keys are chain descriptions and values are lists of letters
        corresponding to that chain description
    letter : string
        letter whose equivalents we are trying to find

    Returns
    -------
    equivalent_letters : list
        list of equivalent chain letters
    """
    equivalent_letters = None
    for k in chain_info:
        for l in chain_info[k]:
            if l == letter:
                equivalent_letters = chain_info[k]
                return equivalent_letters
    print("Failed to find equivalent letters, should not ever be here!")

