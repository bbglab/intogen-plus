import itertools
import logging

from bgreference import refseq


_BUILD = None


def set_build(build):
    """
    Set the build fo the reference genome. Should only be called once

    Args:
        build (str): genome reference build

    """
    global _BUILD
    _BUILD = build
    logging.getLogger(__name__).info('Using %s as reference genome', _BUILD.upper())


def get(chromosome, start, size=1):
    """
    Gets a sequence from the reference genome

    Args:
        chromosome (str): chromosome
        start (int): start position where to look
        size (int): number of bases to retrieve

    Returns:
        str. Sequence from the reference genome

    """
    return refseq(_BUILD, chromosome, start, size)


def get_triplet(chromosome, pos):
    """

    Args:
        chromosome (str): chromosome identifier
        pos (int): central position

    Returns:
        str: 3 bases from the reference genome

    """
    return get(chromosome, pos-1, size=3)


def _slicing_window(seq, n):
    """Yield items of size n through a sequence"""
    it = iter(seq)
    result = ''.join(itertools.islice(it, n))

    if len(result) == n:
        yield result

    for elem in it:
        result = result[1:] + elem
        yield result


def generate_triplets(chromosome, start, stop):
    """Yield triplets from a sequence"""
    seq = get(chromosome, start-1, stop-start+2)
    yield from _slicing_window(seq, 3)
