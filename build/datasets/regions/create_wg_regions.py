import sys

from bgreference import refseq

CHR = [str(i) for i in range(1, 23)] + ['X', 'Y']


def compute_sizes(genome, kmer):
    sizes = []
    for chr_ in CHR:
        seq = refseq(genome, chr_, start=(1 + kmer // 2), size=None)
        sizes.append(tuple(map(str, (chr_, 1 + kmer // 2, len(seq) - kmer // 2))))
    return sizes


def write(sizes):
    print('\t'.join(('CHROMOSOME', 'START', 'END')))
    for s in sizes:
        print('\t'.join(s))


if __name__ == '__main__':
    genome = sys.argv[1]
    kmer = int(sys.argv[2])
    write(compute_sizes(genome, kmer))
