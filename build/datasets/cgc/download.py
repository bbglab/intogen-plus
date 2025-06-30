# Import modules
import gzip
import os
import shutil
import sys
import tarfile

import bglogs
import click
import requests
from homura import download


# Currently version v87
CGC_URL = "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch38/cosmic/v99/Cosmic_CancerGeneCensus_Tsv_v99_GRCh38.tar&bucket=downloads"
TARGET_MEMBER = "Cosmic_CancerGeneCensus_v99_GRCh38.tsv.gz"
COSMIC_KEY = os.getenv("COSMIC_KEY", None)


@click.command()
@click.option('--download', 'download_folder', help='Download folder')
@click.option('--debug', is_flag=True)
def cmdline(download_folder, debug=False):
    bglogs.configure(debug=debug)

    if COSMIC_KEY is None:
        bglogs.error("Environment variable COSMIC_KEY not set")
        bglogs.error("Define your key like this:\n\texport COSMIC_KEY=$(echo \"email@example.com:mycosmicpassword\" | base64)")
        sys.exit(-1)

    output_folder = download_folder
    os.makedirs(output_folder, exist_ok=True)
    output_file = os.path.join(output_folder, 'cancer_gene_census.csv')

    r = requests.get(CGC_URL, headers={"Authorization": "Basic {}".format(COSMIC_KEY)})
    url = r.json()['url']
    bglogs.debug(url)

    if os.path.exists(output_file):
        os.unlink(output_file)

    tar_path = f"{output_file}.tar"
    download(url, path=tar_path)

    # Extract and decompress
    with tarfile.open(tar_path, "r") as tar:
        with tar.extractfile(TARGET_MEMBER) as compressed_file:
            with gzip.open(compressed_file, "rt") as uncompressed:
                with open(output_file, "w") as out_file:
                    shutil.copyfileobj(uncompressed, out_file)


if __name__ == "__main__":
    cmdline()
