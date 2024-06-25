Data Gathering 
---------------

TCGA 
^^^^

TCGA somatic mutations (mc3.v0.2.8 version) were downloaded from
(\ https://gdc.cancer.gov/about-data/publications/pancanatlas\ ). We
then grouped mutations according to their patient’s cancer type into 32
different cohorts. Additionally, we filtered out somatic mutations that
did not pass the somatic filtering from TCGA (i.e., column FILTER !=
“PASS”).

PCAWG
^^^^^

PCAWG somatic mutations were downloaded from the International Cancer
Genome Consortium (ICGC) data portal
(\ https://dcc.icgc.org/releases/PCAWG/consensus_snv_indel/\ ). Note
that only ICGC samples can be freely downloaded from this site, the TCGA
portion of the callsets is controlled data. Instructions on how to
obtain them can be found in the same webpage.

cBioPortal
^^^^^^^^^^

Somatic mutations from Whole Exome Sequencing (WXS) and Whole Genome
Sequencing (WGS) cohorts uploaded in cBioPortal that were not part of
any other projects included in the analysis (i.e., TCGA, PCAWG, ST. JUDE
or HARTWIG) were downloaded (download date 2018/09/01;
http://www.cbioportal.org/datasets). We then created cohorts following
the next criteria:

1. Cohorts with a limited number of samples (i.e., lower than 30
      samples) associated to cancer types with extensive representation
      (such as Breast cancer, Prostate cancer, Colorectal
      adenocarcinoma…) across the compendium of cohorts were removed.

2. Samples were uniquely mapped into one single cohort. If the same
      sample was originally included in two cohorts, we removed the
      sample from one of them.

3. Samples not sequenced from human cancer biopsies were discarded (cell
      lines, xenografts, normal tissue, etc.)

4. When patient information was available, only one sample from the same
      patient was selected. The criteria to prioritize samples from the
      same patient was: WXS over WGS; untreated over treated, primary
      over metastasis or relapse and, finally, by alphabetical order.
      When there is no patient information we assume that all patients
      have only one sample in the cohort.

5. When sequencing platform information was available, samples from the
      same study but with different sequencing platforms were further
      subclassified into WXS and WGS datasets (only if the resulting
      cohorts fulfilled the requirements herein described; otherwise,
      the samples were discarded).

6. When variant calling information was available, samples from the same
      cohort and sequencing type were further classified according to
      their calling algorithm (VarScan, MuTect, etc.). If the resulting
      cohorts for each subclass fulfilled the requirements herein
      described, the samples were included;otherwise, the samples were
      discarded. When variant calling information was not available we
      assumed that all the samples went through the same pipeline.

7. When treatment information was available, samples from the same
      cohort, sequencing type, calling algorithm were further classified
      according to their treatment status (i.e, treated versus
      untreated). If the resulting cohorts from the subclassification
      fulfilled the requirements herein described, the samples were
      included;otherwise, the samples were discarded. When information
      was not available we assumed that samples had not been treated.

8. When biopsy information was available, samples from the same cohort,
      sequencing type, calling algorithm, treatment status were further
      classified according to their biopsy type (i.e, primary, relapse
      or metastasis). If the resulting datasets from the
      subclassification fulfilled the requirements herein described, the
      samples were included; otherwise, the samples were discarded. When
      information was not available we assumed that the biopsy type of
      the sample was primary.

..

HARTWIG
^^^^^^^

Somatic mutations of metastatic WGS from Hartwig Medical Foundation were
downloaded through their platform (download date 2019/05/02). Datasets
were split according to their primary site. Samples from unknown primary
sites (i.e., None, Nan, Unknown, Cup, Na), double primary or aggregating
into cohorts of fewer than 5 samples (i.e., Adrenal, Myeloid, Thymus and
Eye) were not considered. A total of 25 different cohorts were created.

ICGC
^^^^

Somatic mutations from Whole Exome Sequencing (WXS) and Whole Genome
Sequencing (WGS) studies uploaded in ICGC Data Portal
(\ https://dcc.icgc.org/repositories\ ) not overlapping with other
projects included in the analysis (i.e., TCGA, PCAWG, CBIOP or ST. JUDE)
were downloaded on 2018/01/09. We then created cohorts following the
criteria used for the cBioPortal datasets (cBioPortal).

ST. JUDE
^^^^^^^^

Somatic mutations from Whole Exome Sequencing (WXS) and Whole Genome
Sequencing (WGS) of Pediatric Cancer Genome Project uploaded in the St.
Jude Cloud (\ https://www.stjude.cloud/data.html\ ) were downloaded on
2018/07/16. Cohorts were created according to their primary site and
their biopsy type (i.e., primary, metastasis and relapse). Resulting
datasets with fewer than 5 samples were discarded.

PedcBioPortal
^^^^^^^^^^^^^

Somatic mutations from Whole Exome Sequencing (WXS) and Whole Genome
Sequencing (WGS) studies uploaded in PedcBioPortal that were not part of
any other projects included in the analysis (i.e., ST. JUDE or CBIOP)
were downloaded (download date
2018/10/01;\ http://www.pedcbioportal.org/datasets\ ). We then created
cohorts following the criteria described in the cBioPortal dataset
(cBioPortal).

TARGET
^^^^^^

Somatic SNVs from WXS and WGS of two TARGET studies, Neuroblastoma (NB)
and Wilms Tumor (WT), from the TARGET consortium were downloaded from
the GDC Data Portal (download date 2019/03/07).

Beat AML
^^^^^^^^

We downloaded unfiltered somatic mutations from samples included in the
Beat AML study. We next applied the following criteria to create our
Beat AML cohort:

1. We focused on somatic single nucleotide variants from VarScan2 using
      skin as normal control. All samples that did not belong to this
      class were not further analyzed.

2. Samples from relapses were filtered out.

3. Samples from bone-marrow transplants were discarded.

4. If there were several samples per patient fulfilling the points 1-3,
      we selected the first in

..

   chronological order.

257 independent samples of Beat AML tumors composed our Beat AML cohort.

LITERATURE
^^^^^^^^^^

We also manually collected publicly available cohorts from the
literature. Each cohort was filtered following the same steps than
mentioned above for the cBioPortal dataset.

For a summary of all datasets used in the current release of intogen
download Table XX.

.. _section-1:

Preprocessing 
--------------

Given the heterogeneity of the datasets analyzed in the current release
of intOGen (resulting from e.g. differences in the genome aligners,
variant calling algorithms, sequencing coverage, sequencing strategy),
we implemented a robust pre-processing strategy aiming at reducing
possible biases. Specifically, we conducted the following filtering
steps:

1. The pipeline is configured to run using GRCh38 as reference genome.
      Therefore, for each input dataset the pipeline requires that the
      reference genome is defined. Datasets using GRCh37 as reference
      genome were lifted over using PyLiftover
      (\ `https://pypi.org/project/pyliftover <https://pypi.org/project/pyliftover/>`__\ /;
      version 0.3) to GRCh38. Mutations failing to liftover from GRCh37
      to GRCh38 were discarded.

2. We removed mutations with equal alternate and reference alleles,
      duplicated mutations within the sample sample, mutations with ‘N’
      as reference or alternative allele, mutations with a reference
      allele not matching its reference genome and mutations within
      non-canonical chromosomes (i.e., mutations outside chr1 to chr22,
      chrX and chrY).

3. Additionally, we removed mutations with low pileup mappability, i.e.
      mutations in regions that could potentially map elsewhere in the
      genome. For each position of the genome we computed the pileup
      mappability, defined as the average uniqueness of all the possible
      reads of 100bp overlapping a position and allowing up to 2
      mismatches. This value is equal to 1 if all the reads overlapping
      a mutation are uniquely mappable while it is close to 0 if most
      mapping reads can map elsewhere in the genome. Positions with a
      pileup mappability lower than 0.9 were removed from further
      analyses.

4. We filtered out multiple samples from the same donor. The analysis of
      positive selection in tumors requires that each sample in a cohort
      is independent from the other samples. That implies that if the
      input dataset includes multiple samples from the same patient
      --resulting from different biopsy sites, time points or sequencing
      strategies-- the pipeline automatically selects the first
      according to its alphabetical order. Therefore, all mutations in
      the discarded samples are not considered anymore.

5. We also filtered out hypermutated samples. Samples carrying more than
      1000 mutations for WXS and 10000 for WGS and a mutation count
      greater than 1.5 times the interquartile range length above the
      third quartile in their respective dataset were considered
      hypermutated and therefore removed from further analyses.

6. Datasets with filtered synonymous variants are not runnable. Most
      cancer driver identification methods need synonymous variants to
      fit a background mutation model. Therefore, datasets with less
      than 5 synonymous and datasets with a missense/synonymous ratio
      greater than 10 were excluded .

7. When the Variant Effect Predictor (VEP) mapped one mutation into
      multiple transcripts associated with HUGO symbols, we selected the
      MANE transcript of the first HUGO symbol in alphabetical
      order.

8. We also discarded mutations mapping into genes without MANE
      transcript in VEP.111.

Methods for cancer driver gene identification
---------------------------------------------

The current version of the intOGen pipeline uses six cancer driver
identification methods (hereinafter DIMs) to identify cancer driver
genes from somatic point
mutations:\ `dNdScv <https://github.com/im3sanger/dndscv>`__\ and\ `CBaSE <http://genetics.bwh.harvard.edu/cbase/index.html>`__\ ,
which test for mutation count bias in genes against a background
mutation rate model that incorporates regional genomic covariates, the
tumour mutational process and mutation consequence
types;\ `OncodriveCLUSTL <http://bbglab.irbbarcelona.org/oncodriveclustl/home>`__\ ,
which tests for significant clustering of mutations in the protein
sequence; smRegions, which tests for enrichment of mutations in protein
functional domains; HotMAPS, which tests for significant clustering of
mutations in the 3D protein structure;
and\ `OncodriveFML <http://bbglab.irbbarcelona.org/oncodrivefml/home>`__\ ,
which tests for functional impact bias of the observed mutations. Next
we briefly describe the rationale and the configuration used to run each
DIM.

.. _section-2:

|image0|
~~~~~~~~

**dNdScv**

dNdScv [REF] asserts gene-specific positive selection by inferring the
ratio of non-synonymous to synonymous substitutions (dN/dS, hereinafter
termed :math:`\omega`\ ) in the coding region of each gene. The method
resorts to a Poisson-based hierarchical count model that can correct
for: i) the mutational processes operative in the tumor determined by
the mutational profile of single-nucleotide substitutions with its
flanking nucleotides, ii) the distribution of consequence types at the
CDS per gene, and iii) the regional variability of the background
mutation rate; it incorporates information about 10 histone marks from
69 cell lines obtained in ENCODE project [REF].

We downloaded (release date 2023/02/24) and built a new reference
database based on the list of MANE transcripts defined by VEP.111
(GRCh38). We then used this reference database to run dNdScv on all
datasets of somatic mutations using the default setting of the method.

**OncodriveFML**

OncodriveFML [ref] is a tool that aims to detect genes under positive
selection by analysing the functional impact bias of the observed
somatic mutations. Briefly, OncodriveFML consists of three steps: in the
first step, it computes the average Functional Impact (FI) score (in our
pipeline we used CADD v1.4) of coding somatic mutations observed in gene
of interest across a cohort of tumor samples. In the next step, sets of
mutations of the same size as the number of mutations observed in the
gene of interest are randomly sampled following the tri-nucleotide
probabilities. This sampling is repeated N times (N=106 in our
configuration) to generate expected average scores across all genes of
interest. Finally, it compares the observed average FI score with the
expected from the simulations in the form of an empirical p-value. The
p-values are then adjusted with a multiple testing correction using the
Benjamini–Hochberg (FDR).

OncodriveCLUSTL

OncodriveCLUSTL is a sequence-based clustering algorithm to detect
significant linear clustering bias of the observed somatic mutations
[REF]. Briefly, OncodriveCLUSTL first maps somatic single nucleotide
variants observed in a cohort to the genomic element under study. After
smoothing mutations along its genomic sequence using a Tukey kernel
based density function, clusters are detected and scored taking into
account the number and distribution of mutations observed. A score for
each genomic element is obtained by adding up the scores of its
clusters. To estimate the significance of the observed clustering
signals, mutations are locally randomized using tri- or penta-nucleotide
context probabilities calculated from the input cohort.

For this analysis, OncodriveCLUSTL version 1.1.1 was run for the set of
defined MANE transcripts bearing 2 or more SNVs mapping the
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

**CBaSe**

CBaSe [REF] asserts gene-specific positive and negative selection by
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

HotMAPS algorithm (HotMAPS-1.1.3 version REF) was modified to include a
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
(http://karchinlab.org/data/HotMAPS/mupit_modbase.sql.gz). We also
downloaded all protein structures from the PDB (download date
2019/09/20) alongside all human protein 3D models from Modeller
(download date 2019/09/20;
ftp://salilab.org/databases/modbase/projects/genomes/H_sapiens/2013/H_sapiens_2013.tar.xz
and
ftp://salilab.org/databases/modbase/projects/genomes/H_sapiens/2013/ModBase_H_sapiens_2013_refseq.tar.xz).
We then annotated the structures following the steps described in
HotMAPS tutorial
(https://github.com/KarchinLab/HotMAPS/wiki/Tutorial-(Exome-scale)).

Since HotMAPS configuration files are pre-built in GRCh37 coordinates
and our pipeline is designed to run using GRCh38, for each input cohort,
we first converted input somatic mutations to GRCh37, executed the
HotMAPS algorithm and transformed the output coordinates to GRCh38. All
conversions are done using PyLiftover tool.

smRegions

smRegions is a method developed to detect linear enrichment of somatic
mutations in user-defined regions of interest [ref]. Briefly, smRegions
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
positions of all Pfam domains for MANE transcripts (VEP.111) from
BioMart; ii) for each Pfam domain we mapped the first and last amino
acid positions into genomic coordinates using TransVar --using GRCh38 as
reference genome--; iii) we discarded Pfam domains failing to map either
the first or last amino acid positions into genomic coordinates.

smRegions was conceptually inspired by e-driver [ref], although
significant enhancements were introduced. Particularly, i) our
background model accounts for the observed tri-nucleotide frequencies
rather than assuming that all mutations are equally likely; ii) the
statistical test is more conservative; iii) Pfam domains are part of the
required input and can be easily updated by downloading the last Pfam
release iv) the method can be configured to any other setting that aims
to detect genes possibility selected by enrichment of mutations in
pre-defined gene regions.

Combining the outputs of driver identification methods
------------------------------------------------------

Rationale 
~~~~~~~~~~

Our goal is to provide a catalogue of driver elements which
appropriately reflects the consensus from the DIMs we run.

To combine the results of individual statistical tests, p-value
combination methods continue to be a standard approach in the field:
e.g., Fisher [REF], Brown [REF] and Stouffer Z-score [REF] methods have
been used for this purpose. These methods are useful for combining
probabilities in meta-analysis, hence to provide a ranking based on
combined significance under statistical grounds. However, the
application of these methods may bear some caveats:

1. The ranking resulting from p-value combination may lead to
      inconsistencies when compared to the individual rankings, i.e.,
      they may yield a consensus ranking that does not preserve
      recurrent precedence relationships found in the individual
      rankings.

2. Some methods, like Fisher’s or Brown’s method, tend to bear
      anti-conservative performance, thus leading to many likely false
      discoveries.

3. Balanced (non-weighted) p-value combination methods may lead to
      faulty results just because of the influence of one or more DIM
      performing poorly for a given dataset.

Weighted methods to combine p-values, like the weighted Stouffer
Z-score, provide some extra room for proper balancing, in the sense of
modelling the relative credibility of each DIM. We reasoned that any
good operational criteria to allocate weights should satisfy the
following requirements: i) provide a specific weighting for each cohort,
thereby allowing the relative credibility of a DIM to depend on the
cohort; ii) reflect prior knowledge about known bona-fide driver genes;
iii) reflect prior knowledge about the criteria that each DIM employed
to yield its output.

Our approach works independently for each cohort: to create a consensus
list of driver genes for each cohort, we first determine how credible
each DIM is when applied to this specific cohort, on the basis of how
many bona fide cancer genes reported in the COSMIC Cancer Gene Census
database (CGC) are highly ranked according to the DIM. Once the
credibility has been quantified, we use a weighted method for combining
the p-values that each DIM produces for each candidate gene. This
combination takes the DIMs credibility into account. Based on the
combined p-values, we conduct FDR correction to conclude a ranking of
candidate driver genes alongside q-values.

Weight Estimation by Voting
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The relative credibility for each method is based on the ability of the
method to give precedence to well-known genes already collected in the
CGC catalogue of driver genes. As each method yields a ranking of driver
genes, these lists can be combined using a voting system --based on the
so-called Schulze’s method. The method allows us to consider each method
as a voter with some voting rights (weighting) which casts ballots
containing a list of candidates sorted by precedence. Schulze’s method
takes precedence information from each individual method and produces a
new consensus ranking [ref].

Instead of conducting balanced voting, we tune the voting rights of the
methods so that we maximize the enrichment of CGC genes at the top
positions of the consensus list upon voting. In order to prevent
degenerate solutions, we impose some constraints so that every method
contributes with a minimum share. We also limit the combined share of
coalitions of methods based on similar signals of positive selection.
The solution voting rights are deemed the relative credibility for each
method.

Ranking Score
~~~~~~~~~~~~~

Upon selection of a catalogue of bona-fide known driver elements (CGC
catalogue of driver genes) we can provide a score for each ranking
:math:`R` of genes as follows:

:math:`E(R)\  = \frac{p_{i}}{\log(i + 1)}`

where :math:`p_{i}`\ is the proportion of elements with rank higher than
:math:`i` which belong to CGC and N is a suitable threshold to consider
only the N top ranked elements. Using :math:`E` we can define a function
:math:`f` that maps each weighting vector :math:`w` (in the 5-simplex)
to a value :math:`E(R_{w})` where :math:`R_{w}` denotes the consensus
ranking obtained by applying Schulze’s voting with voting rights given
by the weighting vector w.

Optimization with constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Finally we are bound to find a good candidate for
:math:`\widehat{w_{}} = argmax_{}(f)`\ . In order to prevent spurious
outcomes that could be explained by the concordant detection of driver
signals by methods with similar underlying assumptions, we impose a set
of mild ad-hoc constraints. Specifically, we set the following
constraints: i) all relative weights :math:`\geq` 0.05; ii) w(dNdScv) +
w(cBaSe) :math:`\leq` 0.5; iii) w(OncodriveCLUSTL) + w(HotMaps3D) +
w(smRegions) :math:`\leq` 0.5; w(OncodriveFML) :math:`\leq` 0.5.

Optimization is then carried out in two steps. The first step finds a
candidate :math:`\widehat{w_{0}}` by exhaustive search in a rectangular
grid satisfying the constraints defined above (grid step=0.05). In the
second step we take :math:`\widehat{w_{0}}` as the seed for a stochastic
hill-climbing procedure (using Python’s scipy.optimize “basinhopping”,
method=SLSQP and stepsize=0.05).

Estimation of combined p-values using weighted Stouffer Z-score
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using the relative weight estimate that yields a maximum value of the
objective function f we combined the p-values resorting to the weighted
Stouffer Z-score method. Thereafter we performed Benjamini-Hochberg FDR
correction with the resulting combined p-values, yielding one q-value
for each genomic element. If the element belongs to CGC, we computed its
q-value using only the collection of p-values computed for CGC genes
(hereto referred to as CGC q-value). Otherwise, we computed the q-value
using all the computed p-values (hereto referred simply to as q-value).

Tiers of driver genes from sorted list of combined rankings and p-values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To finalize the analysis we considered only genes with at least two
mutated samples in the cohort under analysis. These genes were
classified into four groups according to the level of evidence in that
cohort that the gene harbours positive selection.

1) The first group, named as TIER1, contained genes showing high
      confidence and agreement in their positive selection signals.
      Given the ranked list of genes obtained by the Schulze voting,
      TIER1 comprises all the ranked genes whose ranking is higher than
      the first gene with combined q-value lower than a specific
      threshold (by default threshold=0.05). The second group, name as
      TIER2, was devised to contain known cancer driver genes, showing
      mild signals of positive selection, that were not included in
      TIER1. More in detail, we defined TIER2 genes as those CGC genes,
      not included in TIER2, whose CGC q-value was lower than a given
      threshold (default CGC q-value=0.25). CGC q-value is computed by
      performing multiple test correction of combined p-values
      restricted to CGC genes. The third group, are genes not included
      in TIER1 or TIER2 with scattered signals of positive selection,
      frequently coming from one single method. Particularly, given the
      ranked list of genes by the Schulze voting, TIER3 was composed of
      all the ranked genes with q-value lower than a given threshold (by
      default threshold=0.05) whose ranking is higher than TIER1 last
      gene position and lower than the rejection ranking position. The
      rejection ranking position is defined as the ranking position for
      which all elements have a q-value lower than the input threshold
      (by default threshold=0.05). Finally, other genes not included in
      the aforementioned classes are considered non-driver genes.

|image1|

Combination evaluation
----------------------

We carried out a benchmark of intOGen’s combination against other
procedures in terms of two scoring systems: enrichment of top-ranked
genes in bona-fide known cancer genes (CGC-Score); enrichment in known
false discovery artifacts (Negative-Score), e.g. long, late-replicating,
structural and/or inactive genes.

Datasets for evaluation
~~~~~~~~~~~~~~~~~~~~~~~

We carried out the evaluation across the 32 Whole-Exome cohorts of the
TCGA PanCanAtlas project (downloaded
from\ https://gdc.cancer.gov/about-data/publications/pancanatlas\ ). By
using these cohorts we could ensure that sequence coverage, sequence
alignment and somatic mutation calling was performed using the same
methodology.

The Cancer Genes Census --version v87-- was downloaded from the COSMIC
data portal (\ https://cancer.sanger.ac.uk/census\ ) and used as a
positive set of known cancer driver genes.

We created a catalog of genes that are known not to be involved in
cancerogenesis. This set includes very long genes (HMCN1, TTN, OBSCN,
GPR98, RYR2 and RYR3), and a list of olfactory receptors from Human
Olfactory Receptor Data Exploratorium (HORDE)
(\ https://genome.weizmann.ac.il/horde/\ ). In addition, for all TCGA
cohorts, we added non-expressed genes, defined as genes where at least
80% of the samples showed a negative log2 RSEM. Expression data for TCGA
was downloaded
from\ https://gdc.cancer.gov/about-data/publications/pancanatlas\ .

Scores for evaluation
~~~~~~~~~~~~~~~~~~~~~

We defined a score, so-called CGC-Score, that is intended to measure the
enrichment of CGC elements in the top positions of the ranking;
specifically given a ranking :math:`R` mapping each element to a rank,
we define the CGC-Score of :math:`R` as

:math:`\text{CGC}\text{score}(R)\  = \frac{p_{i}}{log(i + 1)}` /
:math:``\ :math:`\frac{1}{log(i + 1)}`

where :math:`p_{i}`\ is the proportion of elements with rank
:math:`\leq i` that belong to CGC and :math:`N`\ is a suitable threshold
to consider just the top elements in the ranking (by default N=40).

We estimated the CGC-Score across TCGA cohorts for the individual
methods ranking and the combined ranking.

Similarly, we defined another score, so-called Negative-Score, that aims
to measure the proportion of non-cancer genes among the top positions in
the ranking. Particularly, given a ranking :math:`R` mapping each
element to a rank, we define the Negative-Score of :math:`R`\ as

:math:`\text{CGC}\text{score}(R)\  = \frac{p_{i}}{log(i + 1)}` /
:math:``\ :math:`\frac{1}{log(i + 1)}`

where :math:`p_{i}`\ is the proportion of elements with rank
:math:`\leq i` that belong to the negative set and :math:`N` is a
suitable threshold to consider just the top elements in the ranking (by
default N=40). We estimated the Negative-Score across TCGA cohorts for
the individual methods ranking and the combined ranking.

.. _section-3:

Comparison with individual methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When compared to individual methods, the full combination achieves
highest enrichment of CGC genes among top ranked genes. When compared to
commonly used alternative combination methods based on the same
individual outcomes, our proposed combination tended to yield highest
CGC-Score. We also checked the enrichment of known false discovery
artifacts: we found that our combination method tended to yield lower
Negative-Score than individual methods.

Leave-one-out combination
~~~~~~~~~~~~~~~~~~~~~~~~~

We aimed to estimate the contribution of each method’s ranking to the
final consensus ranking. To address this question, we measured the
effect of removing one DIM at a time from the combination by computing
the ratio between the CGC-Score of the combination excluding the DIM and
the GCG-Score of the full combination. This was done for all DIMs across
TCGA cohorts. Only those DIMs that positively contribute to the combined
ranking attained by intOGen show a ratio below one.

We also assessed the effect on the CGC-Score upon leaving each method
out before combining, one method at a time: with few exceptions, the
effect of leaving any method out was detrimental to the CGC-Score,
meaning that in general all the methods are effectively required to
reach the optimal consensus attained. Also, there is substantial
variability across TCGA cohorts in terms of DIM contributions to the
full combination.

In summary, all methods contributed positively to the combinatorial
performance across TCGA supporting our methodological choice for the
individual driver discovery methods.

Comparison with other combinatorial selection methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We compared the CGC-Score and Negative-Score of our combinatorial
selection strategy against other methods frequently used employed to
produce ranking combinations, either based on ranking information --such
as Borda Count (REF)-- or based on statistical information --such as
Fisher (REF) or Brown (REF, REF) methods. Hereto, we briefly describe
the rationale of the four methods we used to benchmark our ranking. For
the sake of compact notation, let’s denote the rank and p-value of gene
:math:`g` produced by method :math:`m_{i}` as :math:`r_{i,\ g}` and
:math:`p_{i,\ g}`\ , respectively.

Borda Count: For each ranked item :math:`g` and method :math:`m_{i},` it
assigns a score :math:`s_{i,\ g} = N - l_{i,\ g},` where :math:`N`
stands for the total number of items to rank and :math:`l_{i,\ g}` is
the number of items ranked below :math:`g` according to method
:math:`m_{i}`\ . For each item :math:`g` an overall score
:math:`s_{g}{\  = \ s}_{1,\ g} + \ldots + s_{k,\ g}` can then be
computed for each :math:`g,` whence a ranking is obtained by descending
sort.

Fisher: It is based on the p-values :math:`p_{i,\ g}`\ . For each item
:math:`g` the method produces a new combined p-value by computing the
statistic:

:math:`F_{g} = - 2\log\ p_{i,\ g}\  \sim \ \chi_{2k}^{2}`\ .

Under the null hypothesis, :math:`F_{g}` are distributed as a chi-square
with :math:`2k` degrees of freedom, whence a p-value, which in turn
yields a raking by ascending sort. Its applicability is limited by the
assumption that the methods provide independent significance tests.

Brown: This method overcomes the independence requirement of Fisher’s
method by modeling the dependencies between the statistical tests
produced by each method, specifically by estimating the covariance
:math:`\Omega_{i,\ j} = cov( - 2\log\text{\ p}_{i,\ g},\  - 2\log\ p_{j,\ g}).`
Brown’s method (REF) and its most recent adaptation (REF) have been
proposed as less biased alternatives to Fisher.

We then computed the CGC-Score and Negative-Score based on the consensus
ranking from the aforementioned combinatorial methods and compared them
to our Schulze’s weighted combination ranking across all TCGA cohorts.
Our combinatorial approach met a larger enrichment in known cancer genes
for 29/32 (90%) TCGA cohorts [Figure].

Drivers postprocessing
----------------------

The intOGen pipeline outputs a ranked list of driver genes for each
input cohort. We aimed to create a comprehensive catalog of driver genes
per tumor type from all the cohorts included in this version.

Then, we performed a filtering on automatically generated driver gene
lists per cohort. This filtering is intended to reduce artifacts from
the cohort-specific driver lists, due to e.g. errors in calling
algorithms, local hypermutation effects, undocumented filtering of
mutations.

We first created a collection of candidate driver genes by selecting
either: i) significant non-CGC genes (q-value < 0.05) with at least two
significant bidders (methods rendering the genes as significant); ii)
significant CGC genes (either q-value < 0.05 or CGC q-value < 0.25) from
individual cohorts. All genes that did not fulfill these requirements
were discarded.

Additionally, candidate driver genes were further filtered using the
following criteria:

1. We discarded non-expressed genes using TCGA expression data. For
      tumor types directly mapping to cohorts from TCGA --including TCGA
      cohorts-- we removed non-expressed genes in that tumor type. We
      used the following criterion for non-expressed genes: genes where
      at least 80% of the samples showed a negative log2 RSEM. For those
      tumor types which could not be mapped to TCGA cohorts this
      filtering step was not done.

2. We also discarded genes highly tolerant to Single Nucleotide
      Polymorphisms (SNP) across human populations. Such genes are more
      susceptible to calling errors and should be taken cautiously. More
      specifically, we downloaded transcript specific constraints from
      gnomAD (release 2.1; 2018/02/14) and used the observed-to-expected
      ratio score (oe) of missense (mys), synonymous (syn) and
      loss-of-function (lof) variants to detect genes highly tolerant to
      SNPs. Genes enriched in SNPs (oe_mys > 1.5 or oe_lof > 1.5 or
      oe_syn > 1.5) with a number of mutations per sample greater than 1
      were discarded. Additionally, we discarded mutations overlapping
      with germline variants (germline count > 5) from a panel of
      normals (PON) from Hartwig Medical Foundation
      (\ https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd/download?path=%2FHMF-Pipeline-Resources&files=SOMATIC_PON.vcf.gz\ ).

3. We also discarded genes that are likely false positives according to
      their known function from the literature. We convened that the
      following genes are likely false positives: i) known long genes
      such as TTN, OBSCN, RYR2, etc.; ii) olfactory receptors from HORDE
      (\ http://bioportal.weizmann.ac.il/HORDE/\ ; download date
      2018/02/14); iii) genes not belonging to Tier1 CGC genes lacking
      literature references according to CancerMine
      (\ http://bionlp.bcgsc.ca/cancermine/\ ).

4. We also removed non CGC genes with more than 3 mutations in one
      sample. This abnormally high number of mutations in a sample may
      be the result of either a local hypermutation process or cross
      contamination from germline variants.

5. Finally we discarded genes whose mutations are likely the result of
      local hypermutation activity. More specifically, some coding
      regions might be the target of mutations associated to COSMIC
      Signature 9 (\ https://cancer.sanger.ac.uk/cosmic/signatures\ )
      which is associated to non-canonical AID activity in lymphoid
      tumours. In those cancer types were Signature 9 is known to play a
      significant mutagenic role (i.e., AML, Non-Hodgkin Lymphomas,
      B-cell Lymphomas, CLL and Myelodysplastic syndromes) we discarded
      genes where more than 50% of mutations in a cohort of patients
      were associated with Signature 9.

Candidate driver genes that were not discarded composed the catalog of
driver genes.

Classification according to annotation level from CGC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We then annotated the catalog of highly confident driver genes according
to their annotation level in CGC. Specifically, we created a three-level
annotation: i) the first level included driver genes with a reported
involvement in the source tumor type according to the CGC; ii) the
second group included CGC genes lacking reported association with the
tumor type; iii) the third group included genes that were not present in
CGC.

To match the tumor type of our analyzed cohorts and the
nomenclature/acronyms of cancer types reported in the CGC we manually
created a dictionary comprising all the names of tumor types from CGC
and cancer types defined in our study, according to the following rules:

1. All the equivalent terms for a cancer type reported in the CGC using
      the Somatic Tumor Type field (e.g. “breast”, “breast carcinoma”,
      “breast cancer”), were mapped into the same tumor type.

2. CGC terms with an unequivocal mapping into our cancer types were
      automatically linked (e.g., “breast” with “BRCA”).

3. CGC terms representing fine tuning classification of a more prevalent
      cancer type that did not represent an independent cohort in our
      study, were mapped to their closest parent tumor type in our study
      (e.g., “malignant melanoma of soft parts” into “cutaneous
      melanoma” or “alveolar soft part sarcoma” into “sarcoma”).

4. Adenomas were mapped to carcinomas of the same cell type
      (e.g.,”hepatic adenoma” into “hepatic adenocarcinoma”, “salivary
      gland adenoma” into “salivary gland adenocarcinoma”).

5. CGC parent terms mapping to several tumor types from our study were
      mapped to each of the potential child tumor types. For instance,
      the term “non small cell lung cancer” was mapped to “LUAD” (lung
      adenocarcinoma) and “LUSC” (lung squamous cell carcinoma).

6. Finally, CGC terms associated with benign lesions, with unspecified
      tumor types (e.g., “other”, “other tumor types”, “other CNS”) or
      with tumor types with missing parent in our study were left
      unmatched.

Mode of action of driver genes 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We computed the mode of action for highly confident driver genes. To do
so, we first performed a pan-cancer run of dNdScv across all TCGA
cohorts. We then applied the aforementioned algorithm (see Mode of
action section below for more details on how the algorithm determines
the role of driver genes according to their distribution of mutations in
a cohort of samples) to classify driver genes into the three possible
roles: Act (activating or oncogene), LoF (loss-of-function or tumor
suppressor) or Amb (ambiguous or non-defined). We then combined these
predictions with prior knowledge from the Cancer Genome Interpreter
[ref] according to the following rules: i) when the inferred mode of
action matched the prior knowledge, we used the consensus mode of
action; ii) when the gene was not included in the prior knowledge list,
we selected the inferred mode of action; iii) when the inferred mode of
action did not match the prior knowledge, we selected that of the prior
knowledge list.

Repository of mutational features
---------------------------------

Linear clusters
~~~~~~~~~~~~~~~

Linear clusters for each gene and cohort were identified by
OncodriveCLUSTL. We defined as significant those clusters in a driver
gene with a p-value lower than 0.05. The start and end of the clusters
were retrieved from the first and last mutated amino acid overlapping
the cluster, respectively.

3D clusters
~~~~~~~~~~~

Information about the positions involved in the 3D clusters defined by
HotMAPS were retrieved from the gene specific output of each cohort. We
defined as significant those amino acids in a driver gene with a q-value
lower than 0.05.

Pfam Domains
~~~~~~~~~~~~

Pfam domains for each driver gene and cohort were identified by
smRegions. We defined as significant those domains in driver genes with
a q-value lower than 0.1 and with positive log ratio of
observed-to-simulated mutations (observed mutations / simulated
mutations > 1). The first and last amino acid are defined from the start
and end of the Pfam domain, respectively.

Excess of mutations
~~~~~~~~~~~~~~~~~~~

The so-called excess of mutations for a given coding consequence-type
quantifies the proportion of observed mutations at this consequence-type
that are not explained by the neutral mutation rate. The excess is
inferred from the dN/dS estimate :math:`\omega` as
:math:`(\omega - 1)\ /\ \omega`\ . We computed the excess for missense,
nonsense and splicing-affecting mutations.

Mode of action
~~~~~~~~~~~~~~

We computed the gene-specific dN/dS estimates for nonsense and missense
mutations, denoted :math:`\omega_{\text{non}}` and
:math:`\omega_{\text{mis}}`\ . Then each gene induces a point in the
plane with coordinates
:math:`(\omega_{\text{non}},\ \omega_{\text{mis}})`\ . We deemed a gene
Act (resp. LoF) if its corresponding point sits above (resp. below) the
diagonal (\ :math:`x = y`\ ) up to an uncertainty threshold (XXX). Genes
within the uncertainty area as well as genes with
:math:`\omega_{\text{non}} < \ 1` and :math:`\omega_{\text{mis}} < \ 1`
were deemed “ambiguous”.

.. |image0| image:: media/image2.png
   :width: 6.5in
   :height: 2.77778in
.. |image1| image:: media/image1.png
   :width: 6.5in
   :height: 1.83333in
