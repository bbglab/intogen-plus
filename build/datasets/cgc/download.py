# Import modules
import os
import sys

import bglogs
import click
import requests
from homura import download


# Currently version v87
CGC_URL = "https://cancer.sanger.ac.uk/cosmic/file_download/GRCh37/cosmic/v87/cancer_gene_census.csv"
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

    download(url, path=output_file)


if __name__ == "__main__":
    cmdline()
