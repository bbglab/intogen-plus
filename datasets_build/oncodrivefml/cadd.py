import gzip
import csv

def count_lines(input_file):
    count = 0
    with gzip.open(input_file, 'rb') as fd:
        for _ in fd:
            count += 1
    return count


def parse_annotations(gene_coordinates, cadd_file):

    total_lines = count_lines(gene_coordinates)

    with gzip.open(gene_coordinates, 'rb') as fd, gzip.open(cadd_file, 'wt') as cadd:
        writer = csv.writer(cadd, delimiter='\t')

        tb = tabix.open(CADD_FILE)
        for line in tqdm(fd, total=total_lines):
            chromosome, status, func, start, end, _, strand, _, info, = line.decode().strip().split('\t')
            if chromosome not in CHROMOSOMES or status != 'protein_coding' or func != 'CDS':
                continue

            for res in tb.querys('{}:{}-{}'.format(chromosome, int(start) - 2, int(end) + 1)):
                _chr, _pos, _ref, _alt, _, _phred = res
                writer.writerow([_chr, _pos, _ref, _alt, _phred])