import gzip
import tarfile
import tempfile
import urllib.request

from pathlib import Path
import click
from tqdm import tqdm


def somatic_filter(input_url, output_file, member_name="SageGermlinePon.98x.38.tsv.gz"):
    temp = tempfile.TemporaryDirectory()

    tmp_tar = Path(temp.name) / "hmf_pipeline_resources.38_v2.0.0--3.tar.gz"
    urllib.request.urlretrieve(input_url, tmp_tar)

    with tarfile.open(tmp_tar, "r:gz") as tar:
        try:
            tar_member_pon = f"hmf_pipeline_resources.38_v2.0--3/dna/variants/{member_name}"
            tmp_pon = Path(temp.name) / tar_member_pon

            tar.extract(tar_member_pon, path=temp.name)
        except KeyError:
            raise FileNotFoundError(f"{member_name} not found in {tmp_tar}")

    with gzip.open(tmp_pon, "rb") as fd, gzip.open(output_file, "wt") as fo:
        for line in tqdm(fd):
            line = line.decode().strip()
            if line.startswith("Chromosome"):
                continue
            chrom, pos, ref, alt, count, _max_reads, _tot_reads = line.split("\t")
            if int(count) > 5:
                fo.write(f"{chrom}\t{pos}\t{ref}\t{alt}\n")

    temp.cleanup()


@click.command()
@click.option("-i", "--input-url", required=True, help="URL of the origina file")
@click.option("-o", "--output-file", type=click.Path(), required=True, help="path of the output file")
def run(input_url, output_file):
    somatic_filter(input_url, output_file)


if __name__ == "__main__":
    run()
