BoostDM Connection
------------------

IntOGen pipeline integrates the generation of files needed to run by BoostDM [1]_ in order to keep a unified data and software environment and limit preprocessing of input for BoostDM as much as possible. The integrations consists in two new steps in the pipeline:

DriverSaturation
^^^^^^^^^^^^^^^^

It computes all the possible mutations for a given gene mapping to the MANE transcript. It uses VEP v111. Specifically, it considers both the exons of the transcript and intronic sites within 25 bps distance from the intron-exon junctions.

Filter MNVs
^^^^^^^^^^^

Individual SNVs in adjacent positions reported in the same sample are discarded as potential multiple nucleotide variants (MNVs) that are wrongly called as separate SNVs.

Parse Profile
^^^^^^^^^^^^^

Profile output by ComputeProfile step is parsed in a 96-channel for every trinucleotide context, to ensure consistency and compatibility with BoostDM run.


.. [1] Ferran Muiños, et al. In silico saturation mutagenesis of cancer genes; Nature 2021. (https://doi.org/10.1038/s41586-021-03771-1)