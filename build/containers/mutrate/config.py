"""Configuration object with all required constants of the project"""


# TODO get rid of this file

import os


class Values:

    """Update the attributes to let the other scripts import them"""

    # dict with site counts for each gene, context and consequence type
    # FIXME missing file. We need to compute it
    SITE_COUNTS_PATH = os.path.join(
        os.environ.get("INTOGEN_DATASETS"), 'shared', 'consequences.pickle.gz'
    )

    # triplet counts genome-wide
    TRI_COUNT_GENOME_PATH = os.path.join(
        os.environ.get("INTOGEN_DATASETS"), 'signature', 'wg.counts.gz'
    )

    # triplet counts exome-wide
    TRI_COUNT_EXOME_PATH = os.path.join(
        os.environ.get("INTOGEN_DATASETS"), 'signature', 'cds.counts.gz'
    )

    # set of signatures
    SIGNATURES_PATH = os.path.join(
        os.environ.get("INTOGEN_DATASETS"), 'deconstructsigs', 'signatures.cosmic.exome.tsv'
    )
