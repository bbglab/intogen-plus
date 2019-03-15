"""
This module contains the methods used to
load and parse the input files: elements and mutations
"""

import gzip
import pickle
import logging
from collections import defaultdict

from bgcache import bgcache
from bgparsers import readers
from intervaltree import IntervalTree

from smregions import reference

logger = logging.getLogger(__name__)


def load_mutations(file):
    """
    Parse the mutations file.
    Mutation in chr M, non SNP, with 'N' in the corresponding triplet and
    with a reference nucleotide different than the one in the reference genome are discarded

    Args:
        file: mutations file
        blacklist (optional): file with blacklisted samples
            Defaults to None.

    Yields:
        One line from the mutations file as a dictionary. Each of the inner elements of
        :ref:`mutations <mutations dict>`

    """
    logger.info('Loading mutations')

    count = 0
    count_discarded = 0
    count_discarded_ref = 0

    for row in readers.variants(file, required=['CHROMOSOME', 'POSITION', 'REF', 'ALT']):
        count += 1
        if row['CHROMOSOME'] == 'M' or row['ALT_TYPE'] != 'snp':
            count_discarded += 1
            continue
        ref_triplet = reference.get_triplet(row['CHROMOSOME'], row['POSITION'])
        if 'N' in ref_triplet or ref_triplet[1] != row['REF']:
            count_discarded_ref += 1
            continue

        yield row

    discarded = round((count_discarded_ref + count_discarded)*100/count, 2)

    if discarded < 25:
        logger.info('Discarded %s %% mutations', discarded)
    elif discarded < 50:
        logger.warning('Discarded %s %% mutations', discarded)
    else:
        logger.warning('Discarded %s %% mutations. Consider revising your mutational dataset or the reference genome you are using.', discarded)


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
    regions_tree = {}
    for i, (k, allr) in enumerate(regions.items()):

        if i % 7332 == 0:
            logger.info("[%d of %d]", i+1, len(regions))

        for r in allr:
            tree = regions_tree.get(r['CHROMOSOME'], IntervalTree())
            tree[r['START']:(r['END']+1)] = r['ELEMENT']
            regions_tree[r['CHROMOSOME']] = tree

    logger.info("[%d of %d]", i+1, len(regions))
    return regions_tree


@bgcache
def load_elements_tree(elements_file):
    elements = readers.elements_dict(elements_file, required=['CHROMOSOME', 'START', 'END', 'ELEMENT'])
    return build_regions_tree(elements)


def load_and_map_variants(variants_file, elements_file, regions_of_interest_file):
    """
    From the elements and variants file, get dictionaries with the segments grouped by element ID and the
    mutations grouped in the same way, as well as some information related to the mutations.

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
    logger.info('Loading elements file')
    elements = readers.elements_dict(elements_file, required=['CHROMOSOME', 'START', 'END', 'ELEMENT', 'SYMBOL'])

    # If the input file is a pickle file do nothing
    if variants_file.endswith(".pickle.gz"):
        with gzip.open(variants_file, 'rb') as fd:
            return pickle.load(fd), elements

    logger.info('Building elements tree')
    elements_tree = load_elements_tree(elements_file)

    # Mapping mutations
    variants_dict = defaultdict(list)
    logger.info("Mapping mutations")
    i = 0
    show_small_progress_at = 100000
    show_big_progress_at = 1000000
    for i, r in enumerate(load_mutations(variants_file)):

        if r['CHROMOSOME'] not in elements_tree:
            continue

        if i % show_small_progress_at == 0:
            print('*', end='', flush=True)

        if i % show_big_progress_at == 0:
            print(' [{} muts]'.format(i), flush=True)

        # Get the interval that include that position in the same chromosome
        intervals = elements_tree[r['CHROMOSOME']][r['POSITION']]

        for interval in intervals:
            element = interval.data
            variants_dict[element].append(r)

    if i > show_small_progress_at:
        print('{} [{} muts]'.format(' '*(((show_big_progress_at-(i % show_big_progress_at)) // show_small_progress_at)+1), i), flush=True)

    regions_of_interest = defaultdict(IntervalTree)
    logger.info("Mapping regions of interest")

    reg_of_interest = readers.elements_dict(regions_of_interest_file)

    for name, data in reg_of_interest.items():
        for r in data:
            if r['CHROMOSOME'] not in elements_tree:
                continue

            intervals = elements_tree[r['CHROMOSOME']][r['START']]

            if len(intervals) == 0:
                logger.warning('Region %s-%s (%s-%s) cannot be mapped' % (r['CHROMOSOME'], r['ELEMENT'], r['START'], r['END']))
            else:
                for interval in intervals:
                    element = interval.data
                    regions_of_interest[element].addi(r['START'], r['END'], r['ELEMENT']+";"+r['SYMBOL'])

    return variants_dict, elements, regions_of_interest
