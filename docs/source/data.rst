Data Gathering
---------------

TCGA
^^^^^

TCGA somatic mutations (mc3.v0.2.8 version) were downloaded from
(`https://gdc.cancer.gov/about-data/publications/pancanatlas <https://gdc.cancer.gov/about-data/publications/pancanatlas>`__).
We then grouped mutations according to their patient’s cancer type into
32 different cohorts. Additionally, we filtered out somatic mutations
that did not pass the somatic filtering from TCGA (i.e., column FILTER
!= “PASS”).

PCAWG
^^^^^

PCAWG somatic mutations were downloaded from the International Cancer
Genome Consortium (ICGC) data portal
(`https://dcc.icgc.org/releases/PCAWG/consensus\_snv\_indel/ <https://dcc.icgc.org/releases/PCAWG/consensus_snv_indel/>`__).
Note that only ICGC samples can be freely downloaded from this site, the
TCGA portion of the callsets is controlled data. Instructions on how to
obtain them can be found in the same webpage.

cBioPortal
^^^^^^^^^^

Somatic mutations from Whole Exome Sequencing (WXS) and Whole Genome
Sequencing (WGS) cohorts uploaded in cBioPortal that were not part of
any other projects included in the analysis (i.e., TCGA, PCAWG, ST. JUDE
or HARTWIG) were downloaded (download date 2018/09/01;
http://www.cbioportal.org/datasets). We then created cohorts following
the next criteria:

1. Cohorts with a limited number of samples (i.e., lower than 30 samples) associated to cancer types with extensive representation (such as Breast cancer, Prostate cancer or Colorectal adenocarcinoma) across the compendium of cohorts were removed.

2. Samples were uniquely mapped into one single cohort. If the same sample was originally included in two cohorts, we removed the sample from one of them.

3. Samples not sequenced from human cancer biopsies were discarded (cell lines, xenografts, normal tissue, etc.).

4. When patient information was available, only one sample from the same patient was selected. The criteria to prioritize samples from the same patient was: WXS over WGS; untreated over treated, primary over metastasis or relapse and, finally, by alphabetical order. When there is no patient information we assume that all patients have only one sample in the cohort.

5. When sequencing platform information was available, samples from the same study but with different sequencing platforms were further subclassified into WXS and WGS datasets (only if the resulting cohorts fulfilled the requirements herein described; otherwise, the samples were discarded).

6. When variant calling information was available, samples from the same cohort and sequencing type were further classified according to their calling algorithm (VarScan, MuTect, etc.). If the resulting cohorts for each subclass fulfilled the requirements herein described, the samples were included;otherwise, the samples were discarded. When variant calling information was not available we assumed that all the samples went through the same pipeline.

7. When treatment information was available, samples from the same cohort, sequencing type, calling algorithm were further classified according to their treatment status (i.e, treated versus untreated). If the resulting cohorts from the subclassification fulfilled the requirements herein described, the samples were included;otherwise, the samples were discarded. When information was not available we assumed that samples had not been treated.

8. When biopsy information was available, samples from the same cohort, sequencing type, calling algorithm, treatment status were further classified according to their biopsy type (i.e, primary, relapse or metastasis). If the resulting datasets from the subclassification fulfilled the requirements herein described, the samples were included; otherwise, the samples were discarded. When information was not available we assumed that the biopsy type of the sample was primary.

Hartwig Medical Foundation
^^^^^^^^^^^^^^^^^^^^^^^^^^

Somatic mutations of metastatic WGS from Hartwig Medical Foundation `https://www.hartwigmedicalfoundation.nl/en/database/ <https://www.hartwigmedicalfoundation.nl/en/database/>`__ were
downloaded through their platform (download date 2019/05/02). Datasets
were split according to their primary site. Samples from unknown primary
sites (i.e., None, Nan, Unknown, Cup, Na), double primary or aggregating
into cohorts of fewer than 5 samples (i.e., Adrenal, Myeloid, Thymus and
Eye) were not considered. A total of 25 different cohorts were created.

ICGC
^^^^

Somatic mutations from Whole Exome Sequencing (WXS) and Whole Genome
Sequencing (WGS) studies uploaded in ICGC Data Portal
(`https://dcc.icgc.org/repositories <https://dcc.icgc.org/repositories>`__)
not overlapping with other projects included in the analysis (i.e.,
TCGA, PCAWG, CBIOP or ST. JUDE) were downloaded on 2018/01/09. We then
created cohorts following the criteria used for the cBioPortal datasets
(cBioPortal).

ST. Jude
^^^^^^^^

Somatic mutations from Whole Exome Sequencing (WXS) and Whole Genome
Sequencing (WGS) of Pediatric Cancer Genome Project uploaded in the St.
Jude Cloud
(`https://www.stjude.cloud/data.html <https://www.stjude.cloud/data.html>`__)
were downloaded on 2018/07/16. Cohorts were created according to their
primary site and their biopsy type (i.e., primary, metastasis and
relapse). Resulting datasets with fewer than 5 samples were discarded.

PedcBioPortal
^^^^^^^^^^^^^

Somatic mutations from Whole Exome Sequencing (WXS) and Whole Genome
Sequencing (WGS) studies uploaded in PedcBioPortal that were not part of
any other projects included in the analysis (i.e., ST. JUDE or CBIOP)
were downloaded (download date 2018/10/01;
`http://www.pedcbioportal.org/datasets <http://www.pedcbioportal.org/datasets>`__).
We then created cohorts following the criteria described in the
cBioPortal dataset (cBioPortal).

TARGET
^^^^^^

Somatic SNVs from WXS and WGS of two TARGET studies, Neuroblastoma (NB)
and Wilms Tumor (WT), from the TARGET consortium were downloaded from
the `Genomic Data Commons Portal <https://gdc.cancer.gov/>`__ (download date 2019/03/07).

Beat AML
^^^^^^^^

We downloaded unfiltered somatic mutations from samples included in the
Beat AML study from the `Genomic Data Commons Portal <https://gdc.cancer.gov/>`__. We next applied the following criteria to create our
Beat AML cohort:

1. We focused on somatic single nucleotide variants from VarScan2 using skin as normal control. All samples that did not belong to this class were not further analyzed.

2. Samples from relapses were filtered out.

3. Samples from bone-marrow transplants were discarded.

4. If there were several samples per patient fulfilling the points 1-3, we selected the first in chronological order.

257 independent samples of Beat AML tumors composed our Beat AML cohort.

Literature
^^^^^^^^^^

We also manually collected publicly available cohorts from the
literature. Each cohort was filtered following the same steps than
mentioned above for the cBioPortal dataset (see above).

.. note:: For further information of all datasets used in the latest release of intOGen, please visit `https://www.intogen.org/beta/download <https://www.intogen.org/beta/download>`__.