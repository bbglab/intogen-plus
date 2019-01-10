"""
This module contains the methods used to
load and parse the input files: elements and mutations

.. _elements dict:

elements (:obj:`dict`)
    contains all the segments related to one element. The information is taken from
    the :file:`elements_file`.
    Basic structure:

    .. code-block:: python

        { element_id:
            [
                {
                'CHROMOSOME': chromosome,
                'START': start_position_of_the_segment,
                'STOP': end_position_of_the_segment,
                'STRAND': strand (+ -> positive | - -> negative)
                'ELEMENT': element_id,
                'SEGMENT': segment_id,
                'SYMBOL': symbol_id
                }
            ]
        }


.. _mutations dict:

mutations (:obj:`dict`)
    contains all the mutations for each element. Most of the information is taken from
    the mutations_file but the *element_id* and the *segment* that are taken from the **elements**.
    More information is added during the execution.
    Basic structure:

    .. code-block:: python

        { element_id:
            [
                {
                'CHROMOSOME': chromosome,
                'POSITION': position_where_the_mutation_occurs,
                'REF': reference_sequence,
                'ALT': alteration_sequence,
                'SAMPLE': sample_id,
                'ALT_TYPE': type_of_the_mutation,
                'CANCER_TYPE': group to which the mutation belongs to,
                'SIGNATURE': a different grouping category,
                }
            ]
        }

.. _mutations data dict:

mutations_data (:obj:`dict`)
    contains the `mutations dict`_ and some metadata information about the mutations.
    Currently, the number of substitutions and indels.
    Basic structure:

    .. code-block:: python

        {
            'data':
                {
                    `mutations dict`_
                },
            'metadata':
                {
                    'snp': amount of SNP mutations
                    'mnp': amount of MNP mutations
                    'mnp_length': total length of the MNP mutations
                    'indel': amount of indels
                }
        }

"""

import gzip
import pickle
import logging
from os.path import exists
from bgcache import bgcache
from bgparsers import readers
from collections import defaultdict
from intervaltree import IntervalTree

from simReg import __logger_name__
from simReg.config import remove_extension_and_replace_special_characters as get_name

logger = logging.getLogger(__logger_name__)


def load_mutations(file, blacklist=None, metadata_dict=None):
    """
    Parsed the mutations file

    Args:
        file: mutations file (see :class:`~simReg.main.OncodriveFML`)
        metadata_dict (dict): dict that the function will fill with useful information
        blacklist (optional): file with blacklisted samples (see :class:`~simReg.main.OncodriveFML`).
            Defaults to None.

    Yields:
        One line from the mutations file as a dictionary. Each of the inner elements of
        :ref:`mutations <mutations dict>`

    """

    # Set of samples to blacklist
    samples_blacklisted = set([s.strip() for s in open(blacklist).readlines()]) if blacklist is not None else set()

    snp = 0

    for row in readers.variants(file, extra=['CANCER_TYPE', 'SIGNATURE'], required=['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE']):
        if row['SAMPLE'] in samples_blacklisted:
            continue
        if row['CHROMOSOME'] == 'M':
            continue

        if row['ALT_TYPE'] == 'snp':
            snp += 1
        else:
            continue

        yield row

    if metadata_dict is not None:
        metadata_dict['snp'] = snp


def build_regions_tree(regions):
    """
    Generates a binary tree with the intervals of the regions

    Args:
        regions (dict): segments grouped by :ref:`elements <elements dict>`.

    Returns:
        dict of :obj:`~intervaltree.IntervalTree`: for each chromosome, it get one :obj:`~intervaltree.IntervalTree` which
        is a binary tree. The leafs are intervals [low limit, high limit) and the value associated with each interval
        is the :obj:`tuple` (element, segment).
        It can be interpreted as:

        .. code-block:: python

            { chromosome:
                (start_position, stop_position +1): (element, segment)
            }

    """
    logger.info('Building regions tree')
    regions_tree = {}
    for i, (k, allr) in enumerate(regions.items()):

        if i % 7332 == 0:
            logger.info("[%d of %d]", i+1, len(regions))

        for r in allr:
            tree = regions_tree.get(r['CHROMOSOME'], IntervalTree())
            tree[r['START']:(r['STOP']+1)] = (r['ELEMENT'], r['SEGMENT'])
            regions_tree[r['CHROMOSOME']] = tree

    logger.info("[%d of %d]", i+1, len(regions))
    return regions_tree


@bgcache
def load_elements_tree(elements_file):
    elements = readers.elements(elements_file)
    return build_regions_tree(elements)


def load_and_map_variants(variants_file, elements_file, regions_of_interest_file, blacklist=None, save_pickle=False):
    """
    From the elements and variants file, get dictionaries with the segments grouped by element ID and the
    mutations grouped in the same way, as well as some information related to the mutations.

    Args:
        variants_file: mutations file (see :class:`~simReg.main.OncodriveFML`)
        elements_file: elements file (see :class:`~simReg.main.OncodriveFML`)
        blacklist (optional): file with blacklisted samples (see :class:`~simReg.main.OncodriveFML`). Defaults to None.
           If the blacklist option is passed, the mutations are not loaded from a pickle file.
        save_pickle (:obj:`bool`, optional): save pickle files

    Returns:
        tuple: mutations and elements

        Elements: `elements dict`_

        Mutations: `mutations data dict`_


    The process is done in 3 steps:
       1. :meth:`load_regions`
       #. :meth:`build_regions_tree`.
       #. each mutation (:meth:`load_mutations`) is associated with the right
          element ID

    """
    # Load elements file
    elements = readers.elements(elements_file)

    # If the input file is a pickle file do nothing
    if variants_file.endswith(".pickle.gz"):
        with gzip.open(variants_file, 'rb') as fd:
            return pickle.load(fd), elements

    # Check if it's already done
    variants_dict_precomputed = variants_file + "_mapping_" + get_name(elements_file) + '.pickle.gz'
    if exists(variants_dict_precomputed) and blacklist is None:
        try:
            logger.info("Using precomputed mutations mapping")
            with gzip.open(variants_dict_precomputed, 'rb') as fd:
                return pickle.load(fd), elements
        except EOFError:
            logger.error("Loading file %s", variants_dict_precomputed)

    # Loading elements tree
    elements_tree = load_elements_tree(elements_file)

    # Mapping mutations
    variants_dict = defaultdict(list)
    variants_metadata_dict = {}
    logger.info("Mapping mutations")
    i = 0
    show_small_progress_at = 100000
    show_big_progress_at = 1000000
    snp_mapped = 0
    for i, r in enumerate(load_mutations(variants_file, metadata_dict=variants_metadata_dict, blacklist=blacklist)):

        if r['CHROMOSOME'] not in elements_tree:
            continue

        if i % show_small_progress_at == 0:
            print('*', end='', flush=True)

        if i % show_big_progress_at == 0:
            print(' [{} muts]'.format(i), flush=True)

        # Get the interval that include that position in the same chromosome
        intervals = elements_tree[r['CHROMOSOME']][r['POSITION']]

        for interval in intervals:
            element, segment = interval.data
            variants_dict[element].append(r)

        if intervals:
            if r['ALT_TYPE'] == 'snp':
                snp_mapped += 1

    if i > show_small_progress_at:
        print('{} [{} muts]'.format(' '*(((show_big_progress_at-(i % show_big_progress_at)) // show_small_progress_at)+1), i), flush=True)

    variants_metadata_dict['snp_mapped'] = snp_mapped
    mutations_data_dict = {'data': variants_dict, 'metadata': variants_metadata_dict}

    if save_pickle:
        # Try to store as precomputed
        try:
            with gzip.open(variants_dict_precomputed, 'wb') as fd:
                pickle.dump(mutations_data_dict, fd)
        except OSError:
            logger.debug("Imposible to write precomputed mutations mapping here: %s", variants_dict_precomputed)

    regions_of_interest = defaultdict(IntervalTree)
    logger.info("Mapping regions of interest")

    reg_of_interest = readers.elements(regions_of_interest_file)

    for name, data in reg_of_interest.items():
        for r in data:
            if r['CHROMOSOME'] not in elements_tree:
                continue

            intervals = elements_tree[r['CHROMOSOME']][r['START']]

            if len(intervals) == 0:
                logger.warning('Region %s-%s (%s-%s) cannot be mapped' % (r['CHROMOSOME'], r['ELEMENT'], r['START'], r['STOP']))
            else:
                for interval in intervals:
                    element, segment = interval.data
                    regions_of_interest[element].addi(r['START'], r['STOP']+1, r['ELEMENT']+";"+r['SYMBOL'])

    logger.info("End regions mapping")

    return mutations_data_dict, elements, regions_of_interest
