import os
import csv
import argparse
import logging
import sqlite3
import time


import pickle
from tqdm import tqdm
from multiprocessing import Pool
from functools import partial
from src import simulate_mutations_signatures, randomizer_aa
from src import utils
from src import simulation_signatures
from src.mutations import mutation_density
from src.pdb_structure import get_structure_info, find_neighbors

logger = logging.getLogger(__name__)  # module logger


def parse_arguments():
    info = 'Detects hotspot protein regions'
    parser = argparse.ArgumentParser(description=info)

    # program arguments
    parser.add_argument('-m', '--mutations',
                        type=str, required=True,
                        help='Mutation counts for specific structures')
    parser.add_argument('-mf', '--maf',
                        type=str, required=True,
                        help='MAF input file')
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
                        default="coordinates.txt.gz",
                        type=str,
                        help='Genomic coordinates of the PDBS')
    parser.add_argument('-sc', '--stop-criterion',
                        default=200,
                        type=int,
                        help='Number of simulations exceeding the maximum observed '
                        'residue before stopping. This speeds computation by spending '
                        'less time on non-significant structures. (Default: 200)')
    parser.add_argument('-c', '--cores',
                        default=4,
                        type=int,
                        help='Maximum processes to run in parallell')
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


def process_structure(structure_id, struct_info, quiet, structure_mutations, df_coordinates, signatures, db):

    # skip structure if no mutations
    if not structure_mutations:
        return None

    # get pdb info
    if 'path' not in struct_info:
        return None

    pdb_path = struct_info.pop('path')

    # read in structure
    structure = utils.read_structure(pdb_path, structure_id, quiet=quiet)
    if structure is None:
        return None

    # make a list of all chain letters in structure
    struct_chains = []
    for k in struct_info.keys():
        struct_chains.extend(struct_info[k])

    # separate out mutation info
    ttypes, mres, mcount, mchains = zip(*structure_mutations)  # if model_mutations else ([], [], [])

    # stratify mutations by their tumor type
    # ttype_ixs is a dictionary that contains
    # ttype as the keys and a list of relevant
    # indices as the values
    unique_ttypes = set(ttypes)
    ttype_ixs = {t: [i for i in range(len(mcount)) if ttypes[i] == t]
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
        return None

    # get neigbours for all residues
    neighbors = find_neighbors(all_res_centers_of_geometry, opts['radius'])
    d_correspondence = simulate_mutations_signatures.generate_correspondence(structure_id, db)
    # iterate through each tumour type
    for tumour in unique_ttypes:
        # skip tumor types if not one specified
        if not opts['tumor_type'] == tumour and not opts['tumor_type'] == 'EVERY':
            continue

        # draw information for the specific tumour type
        t_mut_res_centers_of_geometry = mut_res_centers_of_geometry[tumour]
        t_mut_res_mutation_counts = mut_res_mutation_counts[tumour]

        mut_density = mutation_density(t_mut_res_mutation_counts,
                                                     neighbors)
        mut_vals = mut_density.values()
        if mut_vals:
            max_obs_dens = max(mut_density.values())
        else:
            max_obs_dens = 0

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
        df_protein = df_coordinates[df_coordinates["pdb_id"] == structure_id].copy()
        list_rows = []
        for index, row in df_protein.iterrows():
            list_rows.append(row.to_dict())

        # generate empirical null distribution
        sim_null_dist = simulation_signatures.generate_null_dist_sig(structure_id, list_rows, models, struct_info,
                                                   all_res_centers_of_geometry,
                                                   total_mutations,
                                                   opts['num_simulations'],
                                                   opts['seed'],
                                                   neighbors,
                                                   signatures,
                                                   tumour, d_correspondence,
                                                   opts['stop_criterion'],
                                                   max_obs_dens)
        if len(sim_null_dist) == 0:
            break

        # get a list of lists format for compute p values function
        mut_list = [[res_id, mut_density[res_id]] for res_id in mut_density]
        if not t_mut_res_mutation_counts:
            print("here")

        # aditional information about p-values
        # for specific residues in a structure
        # compute p-values for observed
        obs_pvals, sim_cdf = simulation_signatures.compute_pvals(mut_list, sim_null_dist)

        return [
            structure_id,
            tumour,
            ','.join([str(o[0][1]) for o in mut_list]),
            ','.join([str(o[0][2]) for o in mut_list]),
            ','.join([str(o[0][3][1]) for o in mut_list]),
            ','.join([str(t_mut_res_mutation_counts[o[0]]) for o in mut_list]),
            ','.join([str(o[1]) for o in mut_list]),
            ','.join(map(str, obs_pvals))
        ]


def connect_mysql():

    # make mysql connection
    retries = 5
    while retries > 0:
        try:
            db = sqlite3.connect(os.environ['HOTMAPS_DB'])
            break
        except Exception:
            time.sleep(5)
            retries -= 1

    if retries == 0:
        raise RuntimeError("Impossible to connect")

    return db

def process_structures(quiet, mut_path, file_coordinates, pdb_info):
    output = []
    mutations = utils.read_mutations(mut_path)
    df_coordinates = simulation_signatures.read_file_coordinates(file_coordinates)
    with open(opts['mutations'] + ".signature", "rb") as fd:
        signatures = pickle.load(fd)

    db = connect_mysql()

    for structure_id, struct_info in pdb_info:
        structure_mutations = mutations.get(structure_id, [])

        try:
            result = process_structure(structure_id, struct_info, quiet, structure_mutations, df_coordinates, signatures, db)
        except sqlite3.Error:
            # Try to reconnect DB
            db = connect_mysql()
            result = process_structure(structure_id, struct_info, quiet, structure_mutations, df_coordinates, signatures, db)

        if result is not None:
            output.append(result)

    db.close()

    return output, len(pdb_info)


def main(opts):
    """Currently, performs analysis for the given genes. It attempts to use
    any available PDB sturctures. It then loops through each protein chain
    and tumor type.
    """
    pdb_info = utils.read_pdb_info(opts['annotation'])

    logger.info("Compute signature")
    signatures = randomizer_aa.compute_signature(opts["maf"])
    with open(opts['mutations'] + ".signature", 'wb') as fd:
        pickle.dump(signatures, fd)

    quiet = opts['log_level'] != "DEBUG"

    steps = 100 * opts['cores']
    chunk_size = int(len(pdb_info) / steps) + 1
    process_task = partial(process_structures, quiet, opts['mutations'], opts["genomic_coordinates"])

    header = [[
        'Structure', 'Tumor Type', 'Model', 'Chain', 'Mutation Residues',
        'Residue Mutation Count', 'Mutation Density', 'Hotspot P-value',
    ]]

    with open(opts['output'], 'w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')

        writer.writerows(header)

        # Progress bar
        with tqdm(total=len(pdb_info), desc="Computing PDB structures".rjust(40)) as pb:

            # Multiprocess pool
            pool = Pool(opts['cores'])
            map_method = pool.imap_unordered
            for result, done in map_method(process_task, utils.chunkizator(pdb_info.items(), size=chunk_size)):
                pb.update(done)
                writer.writerows(result)

    logger.info('Done')


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
