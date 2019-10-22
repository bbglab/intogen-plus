"""
Translate COSMIC signatures from genome to CDS context
It takes deconstructSigs format as input: columns ~ A[C>T]G, rows ~ Signature.4
It return the output in the same format.
"""

# Import modules
import os
import sys
import gzip
import json
from shutil import copyfile

import pandas as pd

sys.path.append('../mutrate/')

from utils import deconstruct_to_lex, lex_to_deconstruct, normalize_profile, denorm_profile


INTOGEN_DATASETS = os.path.join(
    '../../datasets', '{}_{}_{}'.format(
        os.environ['INTOGEN_GENOME'],
        os.environ['INTOGEN_VEP'],
        os.environ['INTOGEN_RELEASE']
    )
)
GENOME_TRIPLETS_PATH = os.path.join(INTOGEN_DATASETS, 'shared', 'cds.counts.gz')
EXOME_TRIPLETS_PATH = os.path.join(INTOGEN_DATASETS, 'shared', 'wg.counts.gz')
COSMIC_GENOME_PATH = os.path.join(INTOGEN_DATASETS, 'mutrate', 'signatures.cosmic.genome.tsv')
COSMIC_EXOME_PATH = os.path.join(INTOGEN_DATASETS, 'mutrate', 'signatures.cosmic.exome.tsv')

exists = os.path.isfile(COSMIC_EXOME_PATH)

# Copy the file to the dataset folder
os.makedirs(os.path.dirname(COSMIC_GENOME_PATH), exist_ok=True)
copyfile('signatures.cosmic.genome.tsv', COSMIC_GENOME_PATH)
copyfile('signatures.cosmic.exome.tsv', COSMIC_EXOME_PATH)

# Read the counts
with gzip.open(GENOME_TRIPLETS_PATH, 'rb') as f:
    genome_triplets = json.load(f)

with gzip.open(EXOME_TRIPLETS_PATH, 'rb') as f:
    exome_triplets = json.load(f)

cosmic_genome = pd.read_csv(COSMIC_GENOME_PATH, sep='\t', index_col=0)

# change format of column --context-- labels
cols = list(map(lambda x: deconstruct_to_lex(x), cosmic_genome.columns))
cosmic_genome.columns = cols

data_dict = {c: [] for c in cols}
for sig in cosmic_genome.index:
    profile = dict(zip(cols, cosmic_genome.loc[sig, :].values))
    norm_profile = normalize_profile(profile, genome_triplets)
    norm_profile = denorm_profile(norm_profile, exome_triplets)
    for c in norm_profile:
        data_dict[c].append(norm_profile[c])

df = pd.DataFrame(data_dict, index=cosmic_genome.index)

# change back to decontruct column format
cols = list(map(lambda x: lex_to_deconstruct(x), cols))
df.columns = cols
df.to_csv(COSMIC_EXOME_PATH, sep='\t')
