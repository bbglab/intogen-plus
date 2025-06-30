from collections import namedtuple
import numpy as np
from Bio.PDB import PDBParser
import csv
import itertools as it
import re
import gzip
import src.statistics as mystats


# import modules needed for logging
import logging
import sys
import os
import datetime

logger = logging.getLogger(__name__)  # module logger

# define valid residue names
valid_seqres = set(['ILE', 'GLN', 'GLY', 'GLU', 'CYS',
                    'HIS', 'SER', 'LYS', 'PRO', 'ASN',
                    'VAL', 'THR', 'ASP', 'TRP', 'PHE',
                    'ALA', 'MET', 'LEU', 'ARG', 'TYR',
                    'MSE'])

# pickle header information from read_pickles
header = ["Gene", "mod_base", "RG", "Loops"]
Header = namedtuple('Header', header)
mod_base = ["TargetBeg", "TargetEnd", "SequenceIdentity", "EValue",
            "ga341", "mpqs", "zdope", "PdbCode", "PdbChain",
            "PdbBeg", "PdbEnd", "HitHistory", "tsvmodMethod",
            "tsvmodNo35", "tsvmodRMSD"]
ModBase = namedtuple('ModBase', mod_base)


def start_logging(log_file='', log_level='INFO', verbose=False):
    """Start logging information into the log directory.

    If os.devnull is specified as the log_file then the log file will
    not actually be written to a file.
    """
    if not log_file:
        # create log directory if it doesn't exist
        log_dir = os.path.abspath('log') + '/'
        if not os.path.isdir(log_dir):
            os.mkdir(log_dir)

        # path to new log file
        log_file = log_dir + 'log.run.' + str(datetime.datetime.now()).replace(':', '.') + '.txt'

    # logger options
    lvl = logging.DEBUG if log_level.upper() == 'DEBUG' else logging.INFO
    #myformat = '%(asctime)s - %(name)s - %(levelname)s \n>>>  %(message)s'
    # define logging format
    if verbose:
        myformat = '%(asctime)s - %(name)s - %(levelname)s \n>>>  %(message)s'
    else:
        myformat = '%(message)s'

    # create logger
    if not log_file == 'stdout':
        # normal logging to a regular file
        logging.basicConfig(level=lvl,
                            format=myformat,
                            filename=log_file,
                            filemode='w')
    else:
        # logging to stdout
        root = logging.getLogger()
        root.setLevel(lvl)
        stdout_stream = logging.StreamHandler(sys.stdout)
        stdout_stream.setLevel(lvl)
        formatter = logging.Formatter(myformat)
        stdout_stream.setFormatter(formatter)
        root.addHandler(stdout_stream)
        root.propagate = True


def read_pdb_info(pdb_info_path, pdb_base_path):
    """Reads the pdb info file that relates PDB ID, chain letter, HUGO name,
    and chain description to each other.

    Note: Assumes a multilevel sort of the PDB info file, first by PDB ID and
    second by chain description name.

    Parameters
    ----------
    pdb_info_path : str
        path to PDB info file

    Returns
    -------
    pdb_info_dict : dict
        dict of PDB ID pointing to chain information
    """
    pdb_info_dict = {}
    with open(pdb_info_path) as handle:
        handle.readline()
        myreader = csv.reader(handle, delimiter='\t')
        for pdbid, lines in it.groupby(myreader, lambda x: x[0]):
            # convert from iterator to list
            # for convenience
            lines = list(lines)

            # figure out which chains have the same description
            gene2chain = {}
            for chain_description, lines_subset in it.groupby(lines, lambda x: x[5]):
                gene2chain[chain_description] = [l[1] for l in lines_subset]

            # add path info
            gene2chain['path'] = os.path.join(pdb_base_path, lines[0][4])

            pdb_info_dict[pdbid] = gene2chain
    return pdb_info_dict


def read_mutations(mut_path):
    """Reads the mutation file.

    Note: Assumes a multilevel sort of the PDB info file, first by PDB ID and
    second by chain description name.

    Parameters
    ----------
    mut_path : str
        path to mutation file

    Returns
    -------
    mut_dict : dict
        dict relating structure ID to mutations
    """
    mut_dict = {}
    with open(mut_path) as handle:
        handle.readline()
        myreader = csv.reader(handle, delimiter='\t')
        for pdbid, lines in it.groupby(myreader, lambda x: x[0]):

            # convert lines from an iterator to a list
            # for convenience
            lines = list(lines)

            # currently the regex search is used for odd cases where
            # there is an extra letter with the residue number
            good_ixs = [i for i, m in enumerate(lines)
                        if not re.search('[A-Z]', m[2].split(':')[0])]
            if len(good_ixs) != len(lines):
                logger.debug('Structure {0} has bugged residue numbers'.format(pdbid))
                ####OA PROBLEM HERE ###
                pass

            # keep only residue numbers without problems
            ttype = [lines[i][1] for i in good_ixs]
            mcount = [[int(lines[i][2].split(':')[0]), int(lines[i][3])]
                       for i in good_ixs]
            mutation_chains = [lines[i][2].split(':')[1] for i in good_ixs]  # get chain for pdb

            # reformat output
            mut_info = [[ttype[i], mcount[i][0], mcount[i][1], mutation_chains[i]]
                        for i in range(len(mcount))]

            # save mutations for specific structure
            mut_dict[pdbid] = mut_info
    return mut_dict


def read_structure(pdb_path, structure_id, quiet=True):
    """Reads in a PDB structure.

    Will read gzip compressed PDB structures.

    Parameters
    ----------
    pdb_path : str
        path to pdb file to read
    structure_id : str
        structure id of pdb file

    Returns
    -------
    structure : Bio.PDB structure object | None
       returns PDB structure if possible else none
    """
    pdb_parser = PDBParser(QUIET=quiet)  # parser for pdb files

    # skip if there is no pdb for it
    if not pdb_path:
        logger.debug('Skipping pdb {0}'.format(structure_id))
        return None

    # read in pdb file
    try:
        # handle gziped or uncompressed reading
        if pdb_path.endswith('.gz'):
            with gzip.open(pdb_path, 'rb') as handle:
                structure = pdb_parser.get_structure(structure_id, handle)
        else:
            structure = pdb_parser.get_structure(structure_id, pdb_path)

        # fix homology model chain letters to be "A" instead of " "
        for model in structure:
            for chain in model:
                if chain.id == " ":
                    chain.id = "A"
                    # No need to do this
                    # del model.child_dict[' ']
                    # model.child_dict['A'] = chain

        return structure
    except KeyboardInterrupt:
        # stop if they kill program
        raise
    except:
        logger.info('Fail reading {0}'.format(pdb_path))
        return None


def chunkizator(iterable, size=1000):
    """
    Creates chunks from an iterable

    Args:
        iterable:
        size (int): elements in the chunk

    Returns:
        list. Chunk

    """
    s = 0
    chunk = []
    for i in iterable:
        if s == size:
            yield chunk
            chunk = []
            s = 0
        chunk.append(i)
        s += 1
    yield chunk

