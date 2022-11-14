""" This module can me executed both as a command line and as an imported module

Example of use as a command line:
    python bgscripts/randomizer_aa.py
      -r tests/coordinates_example.txt
      -i tests/mutations_example_filtered.txt
      -o tests/coordinates_example.results
      -s BLCA.signature
      -c BLCA
      -n 10
      --0-based
      --cores=4

Example of use as imported module:
    reg = [
        {'chromosome': '12', 'start': '6122654', 'end': '6122810', 'strand': '-', 'pdb_id': '1atz', 'chain': 'A'},
        {'chromosome': '12', 'start': '6125254', 'end': '6125397', 'strand': '-', 'pdb_id': '1atz', 'chain': 'A'},
        {'chromosome': '12', 'start': '6125681', 'end': '6125821', 'strand': '-', 'pdb_id': '1atz', 'chain': 'A'},
        {'chromosome': '12', 'start': '6125919', 'end': '6126028', 'strand': '-', 'pdb_id': '1atz', 'chain': 'A'}
        ]
    from bgscripts import randomizer_aa
    results = randomizer_aa.randomize_region(
                number_mutations=4, input_regions=reg, number_simulations=1000,
                input_signature='BLCA.signature', cores=4
                )
"""

# Import modules
import argparse
import os
import pickle
import json
import gzip
import csv
from collections import defaultdict, namedtuple
from functools import partial
from itertools import groupby
from concurrent.futures import ProcessPoolExecutor as Pool

import logging
import numpy as np

from bgreference import hg19


# Configure the colorlog module
logger = logging.getLogger()
logger.setLevel(logging.INFO)


# Global variables
Region = namedtuple('Region', 'id chain chromosome start end strand')
# Mutation = namedtuple(
#     'Mutation',
#     ' '.join([
#         'codon', 'alt', 'chromosome', 'position', 'position_upstream',
#         'position_downstream', 'strand', 'aa_ref', 'aa_alt', 'cancer_type',
#         'pdb_id', 'chain'
#         ])
# )
signatures = None
mutations = {}
offset = 0

COMPLEMENTS = dict(zip('ACGTNacgtn', 'TGCANtgcan'))
CODONS = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
    'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
    'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
    'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
    'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
    'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
    'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
    'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
    'GGG': 'G', }
STOP_CODONS = ('TAA', 'TAG', 'TGA')
START_CODONS = ('TTG', 'CTG', 'ATG')


def reverse_complement(sequence):
    """Return the reverse complement of a genomic sequence"""
    seq = list(sequence)
    seq.reverse()
    return ''.join([COMPLEMENTS.get(i, i) for i in seq])


def simulate(items, sample_list, cancer_type=None, simulations=1, cores=1):
    """Read the regions' coordinates and retrieve the corresponding DNA sequence. Then simulates
    missense mutations in accordance with the signature

    :param items: tuple, the first element is the geneid and the second element are the regions
    :param cancer_type: str or None, name of the tumor
    :param simulations: int, number of simulations to do
    :param cores: int, number of cores to use
    :return: tuple, geneid (str) and simulated_mutations (namedtuple)
    """
    geneid, regions = items
    cancer_type = cancer_type if cancer_type is not None else ','.join([str(i) for i in mutations[geneid].keys()])
    changes = []
    prob = dict()
    sequence = ''
    positions = []
    chains = []

    # Group the regions by chain
    sorted_regions = sorted(regions, key=lambda x: x.chain)
    grouped_regions = groupby(sorted_regions, lambda x: x.chain)

    for key, regions in grouped_regions:
        regions = list(regions)
        seq = ''
        pos = []
        chain = []
        strand = set([])
        for i, region in enumerate(sorted(regions, key=lambda x: x.start)):
            start = region.start - 1 if i == 0 else region.start
            end = region.end + 2 if i == len(regions) - 1 else region.end + 1

            if start > end:
                logger.error(
                    'Region {} - {} has a start > end: {} > {}'.format(
                        region.id, region.chain, region.start, region.end
                    )
                )
                continue
            try:
                seq_ = hg19(
                    chromosome=region.chromosome, start=start - offset + 1, size=(end - start)
                )
            except Exception:
                logger.error(
                    'Region {} - {} is outside the chromosome: {}'.format(
                        region.id, region.chain, region
                    )
                )
                continue
            seq += seq_
            for p in range(region.start, region.end + 1):
                pos.append(p)
                chain.append(region.chain)

            strand.add(region.strand)

        # Assumes the strand is the same for each segment of each chain of the region
        if len(strand) > 1:
            logger.error(
                'The chain {} of region {} is annotated with more than one strand'.format(region.chain, region.id)
            )

        if region.strand == '-':
            seq = reverse_complement(seq)
            pos.reverse()

        # Updates the sequence and the list of positions
        sequence += seq
        positions.extend(pos)
        chains.extend(chain)

    codons = [sequence[i: i + 3] for i in range(len(sequence) - 2)]
    for index, (codon, position, chain) in enumerate(zip(codons, positions, chains)):
        # Get rid of codons with 'N'
        if 'N' in codon:
            continue
        alts = ['{}{}{}'.format(codon[0], i, codon[2]) for i in 'ATCG' if i != codon[1]]
        for alt in alts:
            # Only missense variants
            if codon in STOP_CODONS or alt in STOP_CODONS or CODONS[codon] == CODONS[alt]:
                continue
            position_upstream = position - 1 if index == 0 else positions[index - 1]
            position_downstream = position + 1 if index == len(positions) - 1 else positions[index + 1]
            changes.append(
                dict(
                    codon=codon, alt=alt, chromosome=region.chromosome, position=position, strand=region.strand,
                    position_upstream=position_upstream, position_downstream=position_downstream,
                    aa_ref=CODONS[codon], aa_alt=CODONS[alt], cancer_type=cancer_type,
                    pdb_id=geneid, chain=chain,
                )
            )
            ## Edit. Check if type is tuple, otherwise it won't iterate.
            if type(sample_list) is tuple:
                sample_list = list(set(list(sum(sample_list, []))))

            ## Edit. Add sample list to get info for samples.
            for sample in sample_list:
                if sample not in prob:
                    prob[sample] = [signatures[sample]['probabilities'].get((codon, alt), 0) if signatures is not None else 1.0]
                else:
                    prob[sample].append(signatures[sample]['probabilities'].get((codon, alt), 0) if signatures is not None else 1.0)

            #prob.append(signatures['probabilities'].get((codon, alt), 0) if signatures is not None else 1.0)

    # Assumes the length of the sequence is a multiple of 3
    if len(codons) % 3 != 0:
        logger.error(
            'The sequence length of {} is not a multiple of 3: {}'.format(geneid, len(codons))
        )

    # Assumes there are a sufficient number of possible mutations (double of requested)
    num_mutations = sum(mutations[geneid].values())
    if len(changes) < num_mutations * 2:
        logger.warning('There are only {} possible mutations in {}'.format(len(changes), geneid))

    logger.debug('{} - length: {} - missense variants: {}'.format(geneid, len(codons), len(changes) // 3))

    # with Pool(cores) as pool:
    #     fx = partial(
    #         randomize, num_mutations=num_mutations,
    #         p_normalized=p_normalized, changes=changes
    #     )
    #     for simulated_mutations_ in pool.map(fx, range(simulations), chunksize=100):
    #         simulated_mutations.append(simulated_mutations_)
    simulated_mutations = []
    for sample in sample_list:         
        np_prob = np.array(prob[sample])
        p_normalized = np_prob / np.sum(np_prob)
        try:
            simulated_mutations_sample = np.random.choice(
                a=changes,
                size=(simulations, num_mutations),
                p=p_normalized,
                replace=True
            ).flatten().tolist()
        except ValueError:
            return geneid, list()
        simulated_mutations.append(simulated_mutations_sample)

    sim_mut = sum(simulated_mutations, [])
    
    return geneid, sim_mut


# Deprecated!
def randomize(_, num_mutations, p_normalized, changes):
    """Randomize mutations

    :param _: int, index, not used
    :param num_mutations: int, number of mutations to generate
    :param p_normalized: list, normalized signature probabilities
    :param changes: list of namedtuples, possible mutations
    :return:
    """
    simulated_indexes = np.random.choice(range(len(p_normalized)), size=num_mutations, p=p_normalized, replace=True)
    simulated_mutations = {i: x for i, x in enumerate(changes) if i in simulated_indexes}
    simulated_mutations = [simulated_mutations[i] for i in simulated_indexes]
    return simulated_mutations


def load_signature(input_file, load_format):
    """Load precalculated signatures

    :param input_file: path, file with the signature probabilities
    :param load_format: str, pickle or json
    :return: dictionary of signatures
    """
    if os.path.splitext(input_file)[1] == 'gz':
        fx = gzip.open
    else:
        fx = open
    if load_format == 'pickle':
        with fx(input_file, 'rb') as fd:
            signatures = pickle.load(fd)
    elif load_format == 'json':
        with fx(input_file, 'r') as fd:
            signatures = json.load(fd)
    else:
        return None
    return signatures


def compute_signature(mutations_file):

    # Take into account if the mutations are 0 based or 1 based
    offset = 1 #if self.start_at_0 is True else 2
    signatures = {}
    signatures['counts'] = defaultdict(int)

    with open(mutations_file, 'r') as csvfile:
        fd = csv.DictReader(csvfile, delimiter='\t')
        count = 0
        for line in fd:
            chromosome = line["Chromosome"]
            position = int(line["Start_Position"])
            ref = line["Reference_Allele"]
            alt = line["Tumor_Seq_Allele2"]
            if len(ref)!=1 or len(alt)!=1 or ref not in "ACGT" or alt not in "ACGT":
                continue
            signature_ref = hg19(chromosome, position - offset, size=3).upper()
            signature_alt = ''.join([ref[0], alt, ref[-1]])
            signatures['counts'][(signature_ref, signature_alt)] += 1
            count += 1
    signatures['probabilities'] = {k: float(v) / float(count) for k, v in signatures['counts'].items()}

    return signatures


def randomize_region(number_mutations, input_regions, samples, number_simulations=1,
                     start_at_0=True, signature=None, cancer_type=None, cores=1):
    """Randomize a single region instead of a dataset. This modules should be imported in a script.

    :param number_mutations: int, number of mutations to simulate
    :param input_regions: list of tuples, list of lists, or list of dictionaries, with these fields:
           chromosome, start, end, strand, geneid, category
    :param number_simulations: int, number of simulations
    :param input_signature: dict, signatures to use. If None, the mutations will have the same probability to happen
    :param cancer_type: Cancer type, default is None
    :param cores: number of cores to use in the calculation
    :return: list of dictionaries, each dictionary is a simulation

    Example of usage:

    reg = [
        {'chromosome': '12', 'start': '6122654', 'end': '6122810', 'strand': '-', 'pdb_id': '1atz', 'chain': 'A'},
        {'chromosome': '12', 'start': '6125254', 'end': '6125397', 'strand': '-', 'pdb_id': '1atz', 'chain': 'A'},
        {'chromosome': '12', 'start': '6125681', 'end': '6125821', 'strand': '-', 'pdb_id': '1atz', 'chain': 'A'},
        {'chromosome': '12', 'start': '6125919', 'end': '6126028', 'strand': '-', 'pdb_id': '1atz', 'chain': 'A'}
        ]
    from bgscripts import randomizer_aa
    results = randomizer_aa.randomize_region(
                number_mutations=4, input_regions=reg, number_simulations=1000,
                input_signature='BLCA.signature', cores=4
                )

    The output can be read as a pandas dataframe:

    import pandas as pd
    df = pd.DataFrame(results)
    """
    global mutations
    global signatures
    global offset
    offset = 1 if start_at_0 is True else 2

    # Define the regions data
    regions = defaultdict(list)

    if not isinstance(input_regions, list):
        raise TypeError('"input_regions" should be a list of lists, a list of tuples, or a list of dictionaries')

    if len(input_regions) == 0:
        logger.warning('"input_regions" is empty, the region will be skipped')
        return list()

    if isinstance(input_regions[0], dict):
        for reg_ in input_regions:
            regions[reg_['pdb_id']].append(
                Region(
                    id=reg_['pdb_id'], chain=reg_['chain'], chromosome=reg_['chromosome'],
                    start=int(reg_['start']), end=int(reg_['end']), strand=reg_['strand'],
                )
            )
        geneid = reg_['pdb_id']
        # category = reg_['chain']
    elif isinstance(input_regions[0], list) or isinstance(input_regions, tuple):
        for (chromosome, start, end, strand, geneid, chain) in input_regions:
            regions[geneid].append(
                Region(
                    id=geneid, chain=chain, chromosome=chromosome,
                    start=int(start), end=int(end), strand=strand,
                )
            )
    else:
        raise TypeError('"input_regions" should be a list of lists, a list of tuples, or a list of dictionaries')

    # Calculate (or load) the signature
    signatures = signature

    # Prepare the mutation data
    mutations[geneid] = {cancer_type: number_mutations}

    # Simulate the mutations of each region
    results_header = [
        'CHROMOSOME', 'POSITION', 'REF', 'ALT', 'PDB_ID', 'CHAIN', 'STRAND',
        'SIG_REF', 'SIG_ALT', 'AA_REF', 'AA_ALT', 'POSITION_UPSTREAM', 'POSITION_DOWNSTREAM',
        'CANCER_TYPE', 'SIMULATION_ID'
    ]

    results = []
    geneid, simulated_mutations = simulate(
        items=list(regions.items())[0], sample_list=samples, cancer_type=cancer_type, simulations=number_simulations, cores=cores
    )

    for i, mut in enumerate(simulated_mutations):
        ref = mut['codon'][1]
        alt = mut['alt'][1]
        if mut['strand'] == '-':
            ref = COMPLEMENTS[ref]
            alt = COMPLEMENTS[alt]
        results.append(dict(zip(
            results_header, [
                mut['chromosome'], mut['position'], ref, alt, mut['pdb_id'], mut['chain'],
                mut['strand'], mut['codon'], mut['alt'], mut['aa_ref'], mut['aa_alt'],
                mut['position_upstream'], mut['position_downstream'], mut['cancer_type'],
                i % number_simulations
            ]
        )))

    return results


def randomize_dataset(input_mutations, input_regions, input_signature=None, cancer_type=None,
                      output_mutations=None, simulations=1, cores=1, seed=-1):
    """Randomize a whole dataset of regions. This modules is called by the CLI. Writes the results to an output file.

    :param input_mutations: path, file indicating the number of mutations to simulate for each protein. The file has
            the following structure: geneid, cancer_type, number_of_mutations
    :param input_regions: path, file with the coordinates of the regions where the mutations should be generated. The
            file has the following structure: chromosome, start, end, strand, geneid, category
    :param input_signature: path, signatures to use. If None, the mutations will have the same probability to happen
    :param cancer_type: Cancer type, default is None
    :param output_mutations: path, file indicating the number of mutations to simulate for each protein
    :param simulations: int, number of simulations
    :param cores: number of cores to use in the calculation
    :param seed: int, seed to use for the random module
    :return: None
    """
    global mutations
    global signatures

    if output_mutations is None:
        output_mutations = "{}.simulated".format(os.path.basename(input_regions))

    # Set the seed
    if seed >= 0:
        set_seed(seed)

    # Read the regions data
    regions = defaultdict(list)
    with open(input_regions, 'r') as fd:
        for line in fd:
            chromosome, start, end, strand, geneid, category = line.strip().split()
            regions[(geneid, category)].append(
                Region(
                    id=(geneid, category), chromosome=chromosome,
                    start=int(start), end=int(end), strand=strand,
                )
            )

    # Read the mutations data
    if input_mutations is not None:
        muts_ = defaultdict(dict)
        with open(input_mutations, 'r') as fd:
            for line in fd:
                geneid, ct, muts = line.strip().split('\t')
                if cancer_type is None or ct == cancer_type:
                    muts_[geneid][ct] = int(muts)

        for geneid in regions.keys():
            mutations[geneid] = {}
            if geneid[0] in muts_.keys():
                for ct, value in muts_[geneid[0]].items():
                    if ct in mutations[geneid].keys():
                        logger.error('Assigned various number of mutations to the same gene / tumor')
                    mutations[geneid][ct] = value
        del muts_
    else:
        for geneid, values in regions.items():
            length = 0
            for region in values:
                length += region.end - region.start
            mutations[geneid][cancer_type] = length // 100

    # Calculate (or load) the signature
    if input_signature is not None:
        signatures = load_signature(input_signature)

    # Simulate the mutations of each region
    with Pool(cores) as pool:
        with open(output_mutations, 'w', newline='\n') as fd:
            writer = csv.writer(fd, delimiter='\t')
            writer.writerow([
                'CHROMOSOME', 'POSITION', 'REF', 'ALT', 'PROTEIN_ID', 'CHAIN', 'STRAND',
                'SIG_REF', 'SIG_ALT', 'AA_REF', 'AA_ALT', 'POSITION_UPSTREAM', 'POSITION_DOWNSTREAM',
                'CANCER_TYPE', 'SIMULATION_ID'
            ])

            fx = partial(simulate, simulations=simulations, cores=1)
            for geneid, simulated_mutations_ in pool.map(fx, regions.items()):
                for i, simulated_mutations in enumerate(simulated_mutations_):
                    for mut in simulated_mutations:
                        ref = mut.codon[1]
                        alt = mut.alt[1]
                        if mut.strand == '-':
                            ref = COMPLEMENTS[ref]
                            alt = COMPLEMENTS[alt]
                        writer.writerow([
                            mut.chromosome, mut.position, ref, alt, geneid[0], geneid[1],
                            mut.strand, mut.codon, mut.alt, mut.aa_ref, mut.aa_alt,
                            mut.position_upstream, mut.position_downstream, mut.cancer_type, i
                        ])

    return


def set_seed(seed):
    """Set the numpy seed generator

    :param seed: int, seed to use
    :return: None
    """
    np.random.seed(seed)


def cmdline():
    """Parse the command line parser and execute the script

    :return: None
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--regions",
                        dest="input_regions",
                        required=True,
                        help="Input file with region coordinates")
    parser.add_argument("-i", "--input",
                        dest="input_mutations",
                        default=None,
                        help="Input file with the number of mutations that should be generated for each region. " + \
                             "By default it will generate 1 mutation every 100bp.")
    parser.add_argument("-o", "--output",
                        dest="output_mutations",
                        default=None,
                        help="Output file with the randomized mutations")
    parser.add_argument("-c", "--cancer-type",
                        dest="cancer_type",
                        default=None,
                        help="Code of the cancer type")
    parser.add_argument("-s", "--signature",
                        dest="input_signature",
                        default=None,
                        help="Signature to use to simulate mutations")
    parser.add_argument('--0-based',
                        dest='start_at_0',
                        default=False,
                        action='store_true',
                        help="The positions of the mutations are 0 based. Default is False and the positions " +
                             "are considered as 1 based")
    parser.add_argument("-n", "--simulations",
                        dest="simulations",
                        type=int,
                        default=1,
                        help="Number of simulations. Default is 1")
    parser.add_argument("--seed",
                        dest="seed",
                        type=int,
                        default=-1,
                        help="Set the seed generator for the randomization. By default it is chosen randomly")
    parser.add_argument("--cores",
                        dest="cores",
                        type=int,
                        default=os.cpu_count(),
                        choices=list(range(1, os.cpu_count() + 1)),
                        help="Maximum CPU cores to use (default all)")
    args = parser.parse_args()

    global offset
    offset = 1 if args.start_at_0 is True else 2

    randomize_dataset(
        input_mutations=args.input_mutations,    # Number of mutations to generate
        input_regions=args.input_regions,        # Region where to generate mutations
        output_mutations=args.output_mutations,  # File to generate with mutations
        input_signature=args.input_signature,    # Signature to use
        cancer_type=args.cancer_type,
        simulations=args.simulations,            # number of simulation
        cores=args.cores,
        seed=args.seed
    )


if __name__ == "__main__":
    cmdline()
