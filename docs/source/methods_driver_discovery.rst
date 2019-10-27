Methods for cancer driver gene identification
---------------------------------------------

The current version of the intOGen pipeline uses six cancer driver
identification methods (hereinafter DIMs) to identify cancer driver
genes from somatic point mutations:
`dNdScv <https://github.com/im3sanger/dndscv>`__ and
`CBaSE <http://genetics.bwh.harvard.edu/cbase/index.html>`__, which test
for mutation count bias in genes against a background mutation rate
model that incorporates regional genomic covariates, the tumour
mutational process and mutation consequence types;
`OncodriveCLUSTL <http://bbglab.irbbarcelona.org/oncodriveclustl/home>`__,
which tests for significant clustering of mutations in the protein
sequence; smRegions, which tests for enrichment of mutations in protein
functional domains; HotMAPS, which tests for significant clustering of
mutations in the 3D protein structure; and
`OncodriveFML <http://bbglab.irbbarcelona.org/oncodrivefml/home>`__,
which tests for functional impact bias of the observed mutations. Next
we briefly describe the rationale and the configuration used to run each
DIM.


dNdScv
^^^^^^

dNdScv [1]_ asserts gene-specific positive selection by inferring the
ratio of non-synonymous to synonymous substitutions (dN/dS, hereinafter
termed :math:`\omega`) in the coding region of each gene. The method
resorts to a Poisson-based hierarchical count model that can correct
for: i) the mutational processes operative in the tumor determined by
the mutational profile of single-nucleotide substitutions with its
flanking nucleotides, ii) the distribution of consequence types at the
CDS per gene, and iii) the regional variability of the background
mutation rate; it incorporates information about 10 histone marks from
69 cell lines obtained in ENCODE project [2]_.

We downloaded (release date 2018/10/12) and built a new reference
database based on the list canonical transcripts defined by VEP.92
(GRCh38). We then used this reference database to run dNdScv on all
datasets of somatic mutations using the default setting of the method.

OncodriveFML
^^^^^^^^^^^^

OncodriveFML [3]_ is a tool that aims to detect genes under positive
selection by analysing the functional impact bias of the observed
somatic mutations. Briefly, OncodriveFML consists of three steps: in the
first step, it computes the average Functional Impact (FI) score (in our
pipeline we used CADD v1.4) of coding somatic mutations observed in gene
of interest across a cohort of tumor samples. In the next step, sets of
mutations of the same size as the number of mutations observed in the
gene of interest are randomly sampled following the tri-nucleotide
probabilities. This sampling is repeated N times (N=10:sup:`6` in our
configuration) to generate expected average scores across all genes of
interest. Finally, it compares the observed average FI score with the
expected from the simulations in the form of an empirical p-value. The
p-values are then adjusted with a multiple testing correction using the
Benjamini–Hochberg (FDR).

OncodriveCLUSTL
^^^^^^^^^^^^^^

OncodriveCLUSTL is a sequence-based clustering algorithm to detect
significant linear clustering bias of the observed somatic mutations
[4]_. Briefly, OncodriveCLUSTL first maps somatic single nucleotide
variants observed in a cohort to the genomic element under study. After
smoothing mutations along its genomic sequence using a Tukey kernel
based density function, clusters are detected and scored taking into
account the number and distribution of mutations observed. A score for
each genomic element is obtained by adding up the scores of its
clusters. To estimate the significance of the observed clustering
signals, mutations are locally randomized using tri- or penta-nucleotide
context probabilities calculated from the input cohort.

For this analysis, OncodriveCLUSTL version 1.1.1 was run for the set of
defined canonical transcripts bearing 2 or more SNVs mapping the
mutations file as follows: smoothing, clustering windows were kept as
default (11bp). The different consecutive coding sequences contained on
each transcript were concatenated to allow the algorithm to detect
clusters of 2 or more SNVs expanding two exons in a transcript.
Simulations were carried out using previously computed mutational
profiles (see signatures above). All cohorts were run using
tri-nucleotide context SNVs profiles except for cutaneous melanomas,
where penta-nucleotide profiles were calculated. Default randomization
windows of 31bp length were not allowed to expand beyond the coding
sequence boundaries (e.g., windows overlapping part of an exon and an
intron were shifted to fit inside the exon). A total number of 1,000
simulations per transcript were performed. Transcripts with q-value <
0.01 were considered significant.

CBaSe
^^^^^

CBaSe [5]_ asserts gene-specific positive and negative selection by
modelling the distribution of non-synonymous mutation categories under
neutral selection resorting to a Poisson-based hierarchical modelling
approach. As in the case of dNdScv, the method also allows for
correction by i) the mutational processes operative in the tumor as
defined by the mutational profile of single-nucleotide substitutions
--with either tri- or penta- nucleotide context--, ii) the distribution
of consequence types per gene, and iii) regional variability of the
neutral mutation rate. The method theoretically corrects for both known
and unknown covariates of the regional mutation rate, but in practice
the method is sensitive to the count of synonymous mutations. Finally,
the method allows 6 different models based on distinct prior
alternatives for the distribution of the regional mutation rate.

We run a modified version of the CBaSe script to fit the specific needs
of our pipeline. The main modification was adding a clause to
automatically select a regional mutation rate prior distribution that
suits the size of the dataset. Based on the total mutations count in the
dataset, the method runs either an inverse-gamma (mutation count <
12,000), an exponential-inverse-gamma mixture (12,000 < mutation count <
65,000) or a gamma-inverse-gamma mixture (mutation count > 65,000) as
regional mutation rate priors (following communication by Donate
Weghorn, CBaSe’s first author). Furthermore, we also skip the negative
selection analysis and modified the output formatting, including a new
Benjamini-Hochberg FDR column with q-values.

HotMaps3D
^^^^^^^^^

HotMAPS [6]_ algorithm (HotMAPS-1.1.3 version) was modified to include a
new background model that more accurately represents the probability of
somatic mutations in a particular cancer type. The original HotMAPS
algorithm, assumes that all amino-acid substitutions in a protein
structure are equally probable. Herein, we implemented a modified
version of the algorithm where the mutation probability depends on the
trinucleotide context and the mutational processes that are operating in
that cohort of samples. Briefly, for each analyzed protein structure
harbouring missense mutations, the same number of simulated mutations
were randomly generated within the protein structure considering the
precomputed tri-nucleotide frequencies in that cohort. This
pseudo-random sampling was performed N times (N=100,000 in our
configuration) deriving into the background model to compare with the
observed mutational data. The rest of HotMAPS algorithm was kept as it
was.

We downloaded the pre-computed mapping of GRCh37 coordinates into
structure residues from the Protein Data Bank (PDB)
(http://karchinlab.org/data/HotMAPS/mupit\_modbase.sql.gz). We also
downloaded all protein structures from the PDB (download date
2019/09/20) alongside all human protein 3D models from Modeller
(download date 2019/09/20;
ftp://salilab.org/databases/modbase/projects/genomes/H\_sapiens/2013/H\_sapiens\_2013.tar.xz
and
ftp://salilab.org/databases/modbase/projects/genomes/H\_sapiens/2013/ModBase\_H\_sapiens\_2013\_refseq.tar.xz).
We then annotated the structures following the steps described in
HotMAPS tutorial
(https://github.com/KarchinLab/HotMAPS/wiki/Tutorial-(Exome-scale)).

Since HotMAPS configuration files are pre-built in GRCh37 coordinates
and our pipeline is designed to run using GRCh38, for each input cohort,
we first converted input somatic mutations to GRCh37, executed the
HotMAPS algorithm and transformed the output coordinates to GRCh38. All
conversions are done using PyLiftover tool.

smRegions
^^^^^^^^^

smRegions [7]_ is a method developed to detect linear enrichment of somatic
mutations in user-defined regions of interest. Briefly, smRegions
first counts the number of nonsynonymous mutations overlapping with a
Pfam domain in a particular protein. Next, these nonsynonymous variants
are randomly re-sampled N times (N=1,000 in our configuration) along the
nucleotide sequence of the gene, following the probability of mutation
of each base, derived from the pre-computed tri-nucleotide frequencies.
The observed and average number of simulated mutations in the Pfam
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
positions of all Pfam domains for canonical transcripts (VEP.92) from
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

