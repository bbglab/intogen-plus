"""
Translate COSMIC signatures from genome to CDS context
It takes deconstructSigs format as input: columns ~ A[C>T]G, rows ~ Signature.4
It return the output in the same format.
"""

import os
import pandas as pd

import sys
sys.path.append('../mutrate/')

from utils import deconstruct_to_lex, lex_to_deconstruct, normalize_profile, denorm_profile


GENOME_TRIPLETS_PATH = os.path.join(os.environ['INTOGEN_DATASETS'], 'mutrate', 'tri.counts.genome.tsv')
EXOME_TRIPLETS_PATH = os.path.join(os.environ['INTOGEN_DATASETS'], 'mutrate', 'tri.counts.exome.tsv')
COSMIC_GENOME_PATH = os.path.join(os.environ['INTOGEN_DATASETS'], 'mutrate', 'signatures.cosmic.genome.tsv')
COSMIC_EXOME_PATH = os.path.join(os.environ['INTOGEN_DATASETS'], 'mutrate', 'signatures.cosmic.exome.tsv')

exists = os.path.isfile(COSMIC_EXOME_PATH)

genome_triplets = pd.read_csv(GENOME_TRIPLETS_PATH, sep='\t', header=None)
genome_triplets = dict(zip(genome_triplets.iloc[:, 0], genome_triplets.iloc[:, 1]))

exome_triplets = pd.read_csv(EXOME_TRIPLETS_PATH, sep='\t', header=None)
exome_triplets = dict(zip(exome_triplets.iloc[:, 0], exome_triplets.iloc[:, 1]))

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
