
import gzip
import tempfile
import urllib.request

import click
from tqdm import tqdm


def somatic_filter(input_url, output_file):
    temp = tempfile.NamedTemporaryFile()
    urllib.request.urlretrieve(input_url, temp.name)

    with gzip.open(temp.name, 'rb') as fd, gzip.open(output_file, 'wt') as fo:
        for line in tqdm(fd):
            line = line.decode().strip()
            if line.startswith("Chromosome"):
                continue
            chrom, pos, ref, alt, count, _max_reads, _tot_reads = line.split('\t')
            if int(count) > 5:
                fo.write(f'{chrom}\t{pos}\t{ref}\t{alt}\n')

    temp.close()


@click.command()
@click.option('-i', '--input-url', required=True, help='URL of the origina file')
@click.option('-o', '--output-file', type=click.Path(), required=True, help='path of the output file')
def run(input_url, output_file):
    somatic_filter(input_url, output_file)


if __name__ == '__main__':
    run()
