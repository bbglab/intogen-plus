BoostDM Connection
------------------

IntOGen pipeline integrates the generation of files needed to run by BoostDM [1]_ in order to keep a unified data and software environment and limit preprocessing of input for BoostDM as much as possible. The integrations consists in three new steps in the pipeline:

DriverSaturation
^^^^^^^^^^^^^^^^

It computes all the possible mutations for a given gene mapping to the canonical transcript. It uses VEP v101. Specifically, it considers both the exons of the transcript and intronic sites within 25 bps distance from the intron-exon junctions.

Mutrate
^^^^^^^

It scores the mutability of each trinucleotide context (96-channels) from the frequencies observed in the cohort. Since the catalogs of mutations in IntOGen are drawn from possibly different scopes of the genome (e.g. whole genome, whole exome) the frequencies of each trinucleotide context has to be adjusted for the trinucleotide content of the genomic region probed by the sequencing experiment.

Filter MNVs
^^^^^^^^^^^

Individual SNVs in adjacent positions reported in the same sample are discarded as potential multiple nucleotide variants (MNVs) that are wrongly called as separate SNVs.

.. [1] Ferran Muiños, et al. In silico saturation mutagenesis of cancer genes; Nature 2021. (https://doi.org/10.1038/s41586-021-03771-1)