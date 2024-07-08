Methods for cancer driver gene identification
---------------------------------------------

The current version of the intOGen pipeline uses seven cancer driver
identification methods (hereinafter DIMs) to identify cancer driver
genes from somatic point mutations:
`dNdScv <https://github.com/im3sanger/dndscv>`__ and
`cBaSE <http://genetics.bwh.harvard.edu/cbase/index.html>`__, which test
for mutation count bias in genes while correcting for regional genomic
covariates, mutational processes and coding consequence type;
`OncodriveCLUSTL <http://bbglab.irbbarcelona.org/oncodriveclustl/home>`__,
which tests for significant clustering of mutations in the protein
sequence; smRegions, which tests for enrichment of mutations in protein
functional domains; HotMAPS, which tests for significant clustering of
mutations in the 3D protein structure; and
`OncodriveFML <http://bbglab.irbbarcelona.org/oncodrivefml/home>`__,
which tests for functional impact bias of the observed mutations. Next,
we briefly describe the rationale and the configuration used to run each
DIM.


dNdScv
^^^^^^

dNdScv [1]_ asserts gene-specific positive selection by inferring the
ratio of non-synonymous to synonymous substitutions (dN/dS) in the coding
region of each gene. The method resorts to a Poisson-based hierarchical
count model that can correct for: i) the mutational processes operative
in the cohort determined by the mutational profile of single-nucleotide
substitutions with its flanking nucleotides; ii) the regional variability
of the background mutation rate explained by histone modifications -- it
incorporates information about 10 histone marks from 69 cell lines obtained
in ENCODE project [2]_; iii) the abundance of sites per coding consequence
type across in the coding region of each gene.

We downloaded (release date 2023/02/24) and built a new reference
database based on the list of MANE transcripts defined by VEP.111
(GRCh38). We then used this reference database to run dNdScv on all
datasets of somatic mutations using the default setting of the method.

OncodriveFML
^^^^^^^^^^^^

OncodriveFML [3]_ is a tool that aims to detect genes under positive
selection by analysing the functional impact bias of the observed
somatic mutations. Briefly, OncodriveFML consists of three steps: in the
first step, it computes the average Functional Impact (FI) score (in our
pipeline we used CADD v1.6) of coding somatic mutations observed in gene
of interest across a cohort of tumor samples. In the next step, sets of
mutations of the same size as the number of mutations observed in the
gene of interest are randomly sampled following trinucleotide context
conditional probabilities consistent with the relative frequencies of the
mutational profile of the cohort. This sampling is repeated N times
(with N = :math:`10^6` in our configuration) to generate expected average
scores across all genes of interest. Finally, it compares the observed average
FI score with the expected from the simulations in the form of an empirical
p-value. The p-values are then adjusted with a multiple testing correction
using the Benjamini–Hochberg (FDR).

OncodriveCLUSTL
^^^^^^^^^^^^^^^

OncodriveCLUSTL is a sequence-based clustering algorithm to detect
significant linear clustering bias of the observed somatic mutations
[4]_. Briefly, OncodriveCLUSTL first maps somatic single nucleotide
variants (SNVs) observed in a cohort to the genomic element under study. After
smoothing the mutation count per position along its genomic sequence
using a Tukey kernel-based density function, clusters are identified and
scored taking into account the number and distribution of mutations observed.
A score for each genomic element is obtained by adding up the scores of its
clusters. To estimate the significance of the observed clustering
signals, mutations are locally randomized using tri- or penta-nucleotide
context conditional probabilities consistent with the relative frequencies
of the mutational profile of the cohort.

For this analysis, OncodriveCLUSTL version 1.1.3 was run for the set of
defined MANE transcripts bearing 2 or more SNVs mapping the
mutations file. Tuckey-based smoothing was conducted with 11bp windows.
The different consecutive coding sequences contained on
each transcript were concatenated to allow the algorithm to detect
clusters of 2 or more SNVs expanding two exons in a transcript.
Simulations were carried out using pre-computed mutational
profiles. All cohorts were run using tri-nucleotide context SNVs profiles
except for cutaneous melanomas, where penta-nucleotide profiles were calculated.
Default randomization windows of 31bp were not allowed to expand beyond the coding
sequence boundaries (e.g., windows overlapping part of an exon and an
intron were shifted to fit inside the exon). A total number of N = 1,000
simulations per transcript were performed.

cBaSE v1.1
^^^^^

cBaSE [5]_ asserts gene-specific positive and negative selection by
measuring mutation count bias with Poisson-based hierarchical models.
The method allows six different models based on distinct prior
alternatives for the distribution of the regional mutation rate.
As in the case of dNdScv, the method allows for correction by
i) the mutational processes operative in the tumor, with either tri-
or penta- nucleotide context; ii) the site count per consequence type per gene;
iii) regional variability of the neutral mutation rate.

We run a modified version of the cBaSE script to fit the specific needs
of our pipeline. The main modification was adding a rule to automatically
select a regional mutation rate prior distribution. Based on the total
mutation burden in the dataset, the method runs either an inverse-gamma
(mutation count < 12,000), an exponential-inverse-gamma mixture
(12,000 < mutation count < 65,000) or a gamma-inverse-gamma mixture
(mutation count > 65,000) as mutation rate prior distributions -- after
communication with Donate Weghorn, cBaSE’s first author). We also skip the
negative selection analysis part, as it is not needed for downstream analyses.

Mutpanning
^^^^^^^^^^

Mutpanning [9]_ resorts to a mixture signal of positive selection based on two components:
i) the mutational recurrence realized as a Poisson-based count model reminiscent to the
models implemented at dNdScv or cBaSE; ii) a measure of deviance from the characteristic
tri-nucleotide contexts observed in neutral mutagenesis; specifically, an account of the
likelihood that a prescribed count of non-synonymous mutations occur in their observed
given a context-dependent mutational likelihood attributable to the neutral mutagenesis.

HotMaps3D
^^^^^^^^^

HotMAPS [6]_ asserts gene-specific positive selection by measuring
the spatial clustering of mutations in the 3D structure of the protein.
The original HotMAPS method assumes that all amino-acid substitutions in
a protein structure are equally likely. We employed HotMAPS-1.1.3 and
modified it to incorporate a background model that more accurately represents
the mutational processes operative in the cohort.

We implemented a modified version of the method where the trinucleotide
context probability of mutation is compatible with the mutational
processes operative in the cohort. Briefly, for each analyzed protein structure
harbouring missense mutations, the same number of simulated mutations were
randomly generated within the protein structure considering the
precomputed mutation frequencies per tri-nucleotide in the cohort. This
randomization was performed N times (N = :math:`10^5` in our configuration)
thereby leading to a background with which to compare the observed mutational data.
The rest of HotMAPS algorithm was not modified.

We downloaded the pre-computed mapping of GRCh37 coordinates into
structure residues from the Protein Data Bank (PDB)
(`http://karchinlab.org/data/HotMAPS/mupit\_modbase.sql.gz
<http://karchinlab.org/data/HotMAPS/mupit\_modbase.sql.gz>`_).
We also downloaded (on 2019/09/20) all protein structures from the PDB
alongside all human protein 3D models from Modeller
(`ftp://salilab.org/databases/modbase/projects/genomes/H\_sapiens/2013/H\_sapiens\_2013.tar.xz
<ftp://salilab.org/databases/modbase/projects/genomes/H\_sapiens/2013/H\_sapiens\_2013.tar.xz>`_).
and
(`ftp://salilab.org/databases/modbase/projects/genomes/H\_sapiens/2013/ModBase\_H\_sapiens\_2013\_refseq.tar.xz
<ftp://salilab.org/databases/modbase/projects/genomes/H\_sapiens/2013/ModBase\_H\_sapiens\_2013\_refseq.tar.xz>`_).
We then annotated the structures following the steps described in
HotMAPS tutorial (`https://github.com/KarchinLab/HotMAPS/wiki/Tutorial-(Exome-scale)
<https://github.com/KarchinLab/HotMAPS/wiki/Tutorial-(Exome-scale)>`_).

Since HotMAPS configuration files are pre-built in GRCh37 coordinates
and our pipeline is designed to run using GRCh38, for each input cohort,
we first converted input somatic mutations to GRCh37, executed the
HotMAPS algorithm and transformed the output to coordinates to GRCh38. All
conversions were done using the PyLiftover tool.

smRegions
^^^^^^^^^

smRegions [7]_ is a method developed to detect linear enrichment of somatic
mutations in user-defined regions of interest. Briefly, smRegions
first counts the number of non-synonymous mutations overlapping with a
Pfam domain in a particular protein. Next, these non-synonymous variants
are randomized N times (N = 1,000 in our configuration) along the
nucleotide sequence of the gene, following the trinucleotide context
probability derived from precomputed mutation frequencies per tri-nucleotide
in the cohort. The observed and average number of simulated mutations in the Pfam
domain and outside of it are compared using a G-test of goodness-of-fit,
from which the smRegions p-value is derived. We discarded those domains
with a number of observed mutations lower than the average from the
randomizations. The p-values were adjusted with a multiple testing
correction using the Benjamini–Hochberg procedure. Therefore, we
confined the analysis to Pfam domains with a number of observed
mutations higher or equal than the mean simulated number of mutations in
the re-sampling.

To create the database of genomic coordinates of Pfam domains we
followed the next steps: i) we gathered the first and last amino acid
positions of all Pfam domains for MANE transcripts (VEP.111) from
BioMart; ii) for each Pfam domain we mapped the first and last amino
acid positions into genomic coordinates using TransVar --using GRCh38 as
reference genome--; iii) we discarded Pfam domains failing to map either
the first or last amino acid positions into genomic coordinates.

smRegions was conceptually inspired by e-driver [8]_, although
significant enhancements were introduced. Particularly, i) our
background model accounts for the observed tri-nucleotide frequencies
rather than assuming that all mutations are equally likely; ii) the
statistical test is more conservative; iii) Pfam domains are part of the
required input and can be easily updated by downloading the last Pfam
release iv) the method can be configured to any other setting that aims
to detect genes possibility selected by enrichment of mutations in
pre-defined gene regions.


.. [1] Martincorena, I. et al. Universal Patterns of Selection in Cancer and Somatic Tissues. Cell 171, 1029-1041.e21 (2017). doi: 10.1016/j.cell.2017.09.042

.. [2] Roadmap Epigenomics Consortium. Integrative analysis of 111 reference human epigenomes. Nature volume 518, pages 317–330 (19 February 2015). doi: 10.1038/nature14248

.. [3] Loris Mularoni, et al. OncodriveFML: a general framework to identify coding and non-coding regions with cancer driver mutations . Genome Biology (2016)

.. [4] Claudia Arnedo-Pac, et al. OncodriveCLUSTL: a sequence-based clustering method to identify cancer drivers. 2019 Jun 22. Bioinformatics. pii: btz501. doi: 10.1093/bioinformatics/btz501 .

.. [5] Weghorn, et al. D. & Sunyaev, S. Bayesian inference of negative and positive selection in human cancers. Nature Genetics 49, 1785–1788 (2017). doi: 10.1038/ng.3987

.. [6] Tokheim C, et al. Exome-scale discovery of hotspot mutation regions in human cancer using 3D protein structure. Cancer research. 2016a;76:3719–3731. doi: 10.1158/0008-5472.CAN-15-3190

.. [7] Francisco Martínez-Jiménez, et al. Disruption of ubiquitin mediated proteolysis is a widespread mechanism of tumorigenesis. bioRxiv 2019. doi: https://doi.org/10.1101/507764

.. [8] Porta-Pardo E, et al. e-Driver: a novel method to identify protein regions driving cancer. Bioinformatics. 2014;30(21):3109–3114. doi:10.1093/bioinformatics/btu499

.. [9] Dietlein, F., Weghorn, D., Taylor-Weiner, A. et al. Identification of cancer driver genes based on nucleotide context. Nat Genet (2020). https://doi.org/10.1038/s41588-019-0572-y