import argparse
from src.density import *
import src.utils as utils
import numpy as np
import scipy.stats as stats
import csv
import os
import re
from Bio.PDB import *
from src.pdb_structure import *
import src.statistics as mystats
import src.simulation_signatures as sim
import src.simulate_mutations_signatures as gen
import src.simulation as sim2
import src

# import modules needed for logging
import logging
import os

logger = logging.getLogger(__name__)  # module logger

def parse_arguments():
    info = 'Detects hotspot protein regions'
    parser = argparse.ArgumentParser(description=info)

    # program arguments
    parser.add_argument('-m', '--mutations',
                        type=str, required=True,
                        help='Mutation counts for specific structures')
    parser.add_argument('-a', '--annotation',
                        type=str, required=True,
                        help='Annotations about PDB')
    parser.add_argument('-n', '--num-simulations',
                        default=10000,
                        type=int,
                        help='Number of simulations (Default: 10000)')
    parser.add_argument('-r', '--radius',
                        default=10.0,
                        type=float,
                        help='Sphere radius in angstroms (Default: 10.0)')
    parser.add_argument('-s', '--seed',
                        default=101,
                        type=int,
                        help='Random number generator seed (Default: 101)')
    parser.add_argument('-gc', '--genomic_coordinates',
                        default="/home/fran/Documents/clusters3D/clusters/files//coordinates.txt.gz",
                        type=str,
                        help='Genomic coordinates of the PDBS')
    parser.add_argument('-sc', '--stop-criterion',
                        default=200,
                        type=int,
                        help='Number of simulations exceeding the maximum observed '
                        'residue before stopping. This speeds computation by spending '
                        'less time on non-significant structures. (Default: 200)')
    parser.add_argument('-t', '--tumor-type',
                        type=str, default='EVERY',
                        help='Perform analysis for only specific tumor type (Default: "EVERY" = each tumor type)')
    parser.add_argument('-e', '--error-pdb',
                        type=str, default=None,
                        help='File containing structures that have badly formated pdb files')
    
    parser.add_argument('-o', '--output',
                        default='output.txt',
                        type=str,
                        help='Output result file of hotspots')

    # logging arguments
    parser.add_argument('-ll', '--log-level',
                        type=str,
                        action='store',
                        default='',
                        help='Write a log file (--log-level=DEBUG for debug mode, '
                        '--log-level=INFO for info mode)')
    parser.add_argument('-l', '--log',
                        type=str,
                        action='store',
                        default='',
                        help='Path to log file. (accepts "stdout")')
    args = parser.parse_args()

    # handle logging
    if args.log_level or args.log:
        if args.log:
            log_file = args.log
        else:
            log_file = ''  # auto-name the log file
    else:
        log_file = os.devnull
    log_level = args.log_level
    utils.start_logging(log_file=log_file,
                        log_level=log_level)  # start logging

    opts = vars(args)
    return opts


def main(opts):
    """Currently, performs analysis for the given genes. It attempts to use
    any available PDB sturctures. It then loops through each protein chain
    and tumor type.
    """
    # read in data
    logger.info('Reading in annotations . . .')
    pdb_info  = utils.read_pdb_info(opts['annotation'])
    logger.info('Finished reading in annotations.')
    logger.info('Reading in mutations . . .')
    mutations = utils.read_mutations(opts['mutations'])
    logger.info('Finished reading in mutations.')
    df_coordinates = sim.read_file_coordinates(opts["genomic_coordinates"])
    logger.info('Finished reading coordinates')
    # iterate over each structure
    logger.info('Running of PDB structures . . .')
    output = []
    num_pdbs = 0
    num_missing_pdbs = 0
    missing_pdb_list = []
    error_pdb_structs = []
    quiet = True if opts['log_level'] != "DEBUG" else False  # flag indicating pdb warnings
    pdb_parser = PDBParser(QUIET=quiet)  # parser for pdb files
    
    
    for structure_id in pdb_info:
        print (structure_id)
        # get pdb info
        struct_info = pdb_info[structure_id]
        pdb_path = struct_info.pop('path')

        # read in structure
        structure = utils.read_structure(pdb_path, structure_id, quiet=quiet)
        if structure is None:
            continue

        # make a list of all chain letters in structure
        struct_chains = []
        for k in struct_info.keys():
            struct_chains.extend(struct_info[k])

        # get mutation info
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
        unique_ttypes = list(unique_ttypes)

        # obtain relevant info from structure
        tmp_info = get_structure_info(structure, mchains, mres, mcount,
                                      struct_chains, ttype_ixs)
        (mut_res_centers_of_geometry,
         mut_res_mutation_counts,
         all_res_centers_of_geometry,
         models) = tmp_info
        if not all_res_centers_of_geometry:
            logger.error('No available center of geometries for {0}'.format(structure_id))
            continue

        # get neigbours for all residues
        neighbors = find_neighbors(all_res_centers_of_geometry, opts['radius'])
        d_correspondence = gen.generate_correspondence(structure_id)
        # iterate through each tumour type
        for tumour in unique_ttypes:
            # skip tumor types if not one specified
            if (not opts['tumor_type'] == tumour and not opts['tumor_type'] == 'EVERY'):
                continue

            # draw information for the specific tumour type
            t_mut_res_centers_of_geometry = mut_res_centers_of_geometry[tumour]
            t_mut_res_mutation_counts = mut_res_mutation_counts[tumour]

            mut_density = src.mutations.mutation_density(t_mut_res_mutation_counts,
                                                         neighbors)
            mut_vals = mut_density.values()
            if mut_vals:
                max_obs_dens = max(mut_density.values())
            else:
                max_obs_dens =0

            # generate null distribution
            # count total mutations in structure while
            # avoiding double counting due to same id and chain
            # being on multiple models
            obs_models = []
            obs_chains = []
            total_mutations = 0
            for k in t_mut_res_mutation_counts:
                mutations_to_add = t_mut_res_mutation_counts[k]
                for i in range(len(obs_models)):
                    if not k[1] == obs_models[i] and k[2] == obs_chains[i]:
                        mutations_to_add = 0
                        break
                total_mutations += mutations_to_add
                obs_models.append(k[1])
                obs_chains.append(k[2])
            df_protein = df_coordinates[df_coordinates["pdb_id"]==structure_id].copy()
            list_rows = []
            for index,row in df_protein.iterrows():
                list_rows.append(row.to_dict())
            
            # generate empirical null distribution
            print total_mutations
            sim_null_dist = sim.generate_null_dist_sig(structure_id,list_rows,models,struct_info,                                              all_res_centers_of_geometry,
                                                   total_mutations,
                                                   opts['num_simulations'],
                                                   opts['seed'],
                                                   neighbors,
                                                   os.environ['DATASETS_FOLDER']+"/signatures/"+tumour+".signature",
                                                   tumour,d_correspondence,
                                                   opts['stop_criterion'],
                                                   max_obs_dens)
            if len(sim_null_dist) == 0:
            	break
	    '''
            sim_null_dist = sim2.generate_null_dist(structure_id, models, struct_info,
                                                   all_res_centers_of_geometry,
                                                   total_mutations,
                                                   opts['num_simulations'],
                                                   opts['seed'],
                                                   neighbors,
                                                   opts['stop_criterion'],
                                                   max_obs_dens)
                                                           
            '''
            # get a list of lists format for compute p values function
            mut_list = [[res_id, mut_density[res_id]] for res_id in mut_density]
            if not t_mut_res_mutation_counts:
                print("here")

            # aditional information about p-values
            # for specific residues in a structure
            # compute p-values for observed
            obs_pvals, sim_cdf = sim.compute_pvals(mut_list, sim_null_dist)

            output.append([structure_id, tumour,
                            ','.join([str(o[0][1]) for o in mut_list]),
                            ','.join([str(o[0][2]) for o in mut_list]),
                            ','.join([str(o[0][3][1]) for o in mut_list]),
                            ','.join([str(t_mut_res_mutation_counts[o[0]])
                                        for o in mut_list]),
                            ','.join([str(o[1]) for o in mut_list]),
                            ','.join(map(str, obs_pvals)),])

    # write output to file
    output = [['Structure', 'Tumor Type', 'Model', 'Chain', 'Mutation Residues',
               'Residue Mutation Count', 'Mutation Density', 'Hotspot P-value',
              ]] + output
    with open(opts['output'], 'w') as handle:
        csv.writer(handle, delimiter='\t', lineterminator='\n').writerows(output)

    # if user specified to log failed reading of pdbs
    if opts['error_pdb'] and error_pdb_structs:
        with open(opts['error_pdb'], 'w') as handle:
            for bad_pdb in error_pdb_structs:
                handle.write(bad_pdb+'\n')

    print("NUM_MODEL_DIFF: " + str(sim.NUM_MODEL_DIFF))
    print("NUM_CHAIN_DIFF: " + str(sim.NUM_CHAIN_DIFF))
    print("STRUCT_MODEL_DIFF: " + str(sim.STRUCT_MODEL_DIFF))
    print("STRUCT_CHAIN_DIFF: " + str(sim.STRUCT_CHAIN_DIFF))
    logger.info('Finished successfully!')


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
