from src import utils
from src import pdb_structure as pstruct

# get logger
import logging
import os
logger = logging.getLogger(__name__)  # module logger


def mutation_density(mut_ct_dict, neighbours):
    """
    Find the mutation density at each residue position

    Parameters
    ----------
    mut_ct_dict: dict
        dictionary containing full residue ids as keys
        and mutation counts as values
    neighbours: dict
        dictionary containing full residue ids as keys
        and lists of neigbours' residue ids as values

    Returns
    -------
    density_dict : dict
        Maps each residue to a number of mutations within
        its neighbors and itself. Neighbors are defined
        by an angstrom radius cutoff.
    """

    density_dict = {}
    for Residue1 in mut_ct_dict:
        # count muations on the residue itself
        num_mutations = mut_ct_dict[Residue1]
        # count mutations on neighbouring residues
        for Residue2 in neighbours[Residue1]:
            # will be present in mut_ct_dict
            # only if there was a mutation
            if Residue2 in mut_ct_dict:
                num_mutations += mut_ct_dict[Residue2]
        density_dict[Residue1] = num_mutations

    return density_dict


def summarize_residues(mutations, pdb_info, radius,
                       rASA, dssp, tmp_dir,
                       quiet=True):
    # iterate over each structure
    logger.info('Running of PDB structures . . .')
    output = [['structure', 'tumor type', '# buried residues',
               '# protein interface residues', '# nucleic acid interface residues',
               'total residues', '# buried mutations', '# protien interface mutations',
               '# nucleic acid interface mutations', 'total # mutations',
               'burial p-value', 'protein interface p-value',
               'nucleic acid interface p-value']]
    for structure_id in pdb_info:
        #if structure_id.startswith('ENSP') or structure_id.startswith('NP_'):
            #continue
        #print structure_id
        # get pdb info
        struct_info = pdb_info[structure_id]

        if 'path' not in struct_info:
            continue

        pdb_path = struct_info.pop('path')

        # read in structure
        structure = utils.read_structure(pdb_path, structure_id, quiet=quiet)
        if structure is None:
            continue

        # make a list of all chain letters in structure
        struct_chains = []
        for k in struct_info.keys():
            struct_chains.extend(struct_info[k])

        structure_mutations = mutations.get(structure_id, [])
        # skip structure if no mutations
        if not structure_mutations:
            continue

        # separate out mutation info
        ttypes, mres, mcount, mchains = zip(*structure_mutations) # if model_mutations else ([], [], [])

        # stratify mutations by their tumor type
        # ttype_ixs is a dictionary that contains
        # ttype as the keys and a list of relevant
        # indices as the values
        unique_ttypes = set(ttypes)
        ttype_ixs = {t: [i for i in range(len(mcount)) if ttypes[i]==t]
                     for t in unique_ttypes}
        #ttype_ixs['PANCAN'] = range(len(mcount))
        # add PANCAN as a "tumour type"
        unique_ttypes = list(unique_ttypes)
        #unique_ttypes.append('PANCAN')

        # obtain relevant info from structure
        tmp_info = pstruct.get_structure_info(structure, mchains, mres, mcount,
                                              struct_chains, ttype_ixs)
        (mut_res_centers_of_geometry,
         mut_res_mutation_counts,
         all_res_centers_of_geometry,
         models) = tmp_info

        annotated_chains = {chain
                            for description in struct_info
                            for chain in struct_info[description]}

        # find buried residues
        buried_res = pstruct.get_buried_residues(structure, rASA, tmp_dir, dssp)
        tmp_buried = [res_id
                      for res_id in buried_res
                      if res_id[2] in annotated_chains]
        total_res = len(tmp_buried)
        buried_res_info = {(info[1], info[2], info[3])
                           for info in tmp_buried
                           if info[-1] == 1}
        num_buried_res = len(buried_res_info)

        # find interface residues for proteins and nucleic acids
        interface_res = pstruct.get_interface_residues(structure, radius)
        interface_prot_info = {(res_id[1], res_id[2], res_id[3][1])
                               for res_id in interface_res
                               if (res_id[2] in annotated_chains) and interface_res[res_id][0]==1}
        interface_na_info = {(res_id[1], res_id[2], res_id[3][1])
                             for res_id in interface_res
                             if (res_id[2] in annotated_chains) and sum(interface_res[res_id][1:])>=1}
        num_interface_prot_res = len(interface_prot_info)
        num_interface_na_res = len(interface_na_info)

        # iterate through each tumour type
        pan_counts = []
        pan_buried_counts = []
        pan_interface_prot_counts, pan_interface_na_counts = [], []
        tmp_output = []
        for tumour in unique_ttypes:
            # skip tumor types if not one specified
            #if (not opts['tumor_type'] == tumour and not opts['tumor_type'] == 'EVERY'):
                #continue

            # draw information for the specific tumour type
            t_mut_res_centers_of_geometry = mut_res_centers_of_geometry[tumour]
            t_mut_res_mutation_counts = mut_res_mutation_counts[tumour]

            # count total mutations in structure while
            # avoiding double counting due to same id and chain
            # being on multiple models
            obs_models = []
            obs_chains = []
            total_mutations = 0
            total_buried_muts = 0
            total_interface_prot_muts, total_interface_na_muts = 0, 0
            banned_chains = set()
            #if not tumour == 'PANCAN':
            if True:
                for k in t_mut_res_mutation_counts:
                    mutations_to_add = t_mut_res_mutation_counts[k]

                    # prevent double counting
                    cur_model = k[1]
                    cur_chain = k[2]
                    cur_pos = k[3][1]
                    for i in range(len(obs_models)):
                        if not cur_model == obs_models[i] and cur_chain == obs_chains[i]:
                            mutations_to_add = 0
                            break
                    if (cur_chain, cur_pos) in banned_chains:
                        mutations_to_add = 0

                    # add all equivalent chains to banned list
                    equiv_chains = pstruct.find_eq_letters(struct_info, cur_chain)
                    if equiv_chains is not None:
                        equiv_pos = set([(e, cur_pos) for e in equiv_chains])
                        banned_chains |= equiv_pos - set([(cur_chain, cur_pos)])

                    # add to total mutation count
                    total_mutations += mutations_to_add

                    # current residue of interest
                    curr_res = (cur_model, cur_chain, cur_pos)

                    # add buried residue mutation counts
                    is_buried = [(m, c[0], c[1]) in buried_res_info
                                 for c in equiv_pos
                                 for m in range(4)]
                    #if (curr_res in buried_res_info):
                    if any(is_buried):
                        total_buried_muts += mutations_to_add
                        pan_buried_counts.append(mutations_to_add)

                    # add interface residue mutation counts
                    is_interface_prot = [(m, c[0], c[1]) in interface_prot_info
                                         for c in equiv_pos
                                         for m in range(4)]
                    is_interface_na = [(m, c[0], c[1]) in interface_na_info
                                       for c in equiv_pos
                                       for m in range(4)]
                    #if (curr_res in interface_info):
                    if any(is_interface_prot):
                        total_interface_prot_muts += mutations_to_add
                        pan_interface_prot_counts.append(mutations_to_add)
                    if any(is_interface_na):
                        total_interface_na_muts += mutations_to_add
                        pan_interface_na_counts.append(mutations_to_add)

                    # mark chains/models
                    obs_models.append(k[1])
                    obs_chains.append(k[2])
                pan_counts.append(total_mutations)
            else:
                total_mutations = sum(pan_counts)
                total_buried_muts = sum(pan_buried_counts)
                total_interface_prot_muts = sum(pan_interface_prot_counts)
                total_interface_na_muts = sum(pan_interface_na_counts)

            tmp_output.append([structure_id, tumour, num_buried_res,
                               num_interface_prot_res, num_interface_na_res,
                               total_res, total_buried_muts,
                               total_interface_prot_muts, total_interface_na_muts,
                               total_mutations,
                               ])

        output.extend(tmp_output)
    return output
