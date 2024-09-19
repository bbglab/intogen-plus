Data Collection
---------------

TCGA
^^^^^

TCGA somatic mutations (mc3.v0.2.8 version) were downloaded from
(`https://gdc.cancer.gov/about-data/publications/pancanatlas <https://gdc.cancer.gov/about-data/publications/pancanatlas>`__).
We then grouped mutations according to their patient’s cancer type into
32 different cohorts. Additionally, we kept somatic mutations
passing the somatic filtering from TCGA (i.e., column FILTER
== “PASS”).

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
any other projects included in the analysis (i.e., TCGA, PCAWG, St. Jude
or HARTWIG) were downloaded on 2018/09/01 (http://www.cbioportal.org/datasets).
We then created cohorts following the next criteria:

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
downloaded on 2023/09/14 through their platform. Datasets
were split according to their DOID identification and linked to an oncotree identification. Samples from unknown primary
sites (i.e., None, Nan, Unknown, Cup, Na) or  double primary  were not considered. When a tumor type groups less than 20 samples, parental oncotree levels were considered up until grouping at least 20 samples per cohort. A total of 38 different cohorts were created. 

ICGC
^^^^

Somatic mutations from Whole Exome Sequencing (WXS) and Whole Genome
Sequencing (WGS) studies uploaded in ICGC Data Portal
(`https://dcc.icgc.org/repositories <https://dcc.icgc.org/repositories>`__)
not overlapping with other projects included in the analysis (i.e.,
TCGA, PCAWG, CBIOP or St. Jude) were downloaded from release 2019 on on 2021/10/21. We then
created cohorts following the criteria used for the cBioPortal datasets
(cBioPortal).

St. Jude
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
any other projects included in the analysis (i.e., St. Jude or CBIOP)
were downloaded on 2018/10/01 (`http://www.pedcbioportal.org/datasets <http://www.pedcbioportal.org/datasets>`__).
We then created cohorts following the criteria described in the
cBioPortal dataset (cBioPortal).

TARGET
^^^^^^

Somatic SNVs from WXS and WGS of four TARGET studies, Neuroblastoma (NB)
and Wilms Tumor (WT), from the TARGET consortium were downloaded on 2019/03/07, 
Osteosarcoma (OS) and Acute Myeloid Leukemia (AML) we downloaded in 2020/09/17  from
the `Genomic Data Commons Portal <https://gdc.cancer.gov/>`__.

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

CGCI
^^^^
Somatic mutations from Whole Genome Sequencing (WGS) of the The Cancer Genome Characterization Initiative (`CGCI 
<https://www.cancer.gov/ccg/research/genome-sequencing/cgci>`__) were downloaded from the Genomic Data Commons (`GDC 
portal <https://portal.gdc.cancer.gov/projects>`__) on 2021/05/06.

CPTAC
^^^^^
Somatic mutations from Whole Exome Sequencing (WXS)  of the Clinical Proteomic Tumor Analysis Consortium (`CPTAC 
<https://proteomics.cancer.gov/programs/cptac>`__) were downloaded from the `GDC portal <https://portal.gdc.cancer.gov/projects>`__ on 2021/05/06.

Literature
^^^^^^^^^^

We also manually collected publicly available cohorts from the
literature. Each cohort was filtered following the same steps than
mentioned above for the cBioPortal dataset (see above).

.. note:: For further information of all datasets used in the latest release of intOGen, please visit `https://www.intogen.org/beta/download <https://www.intogen.org/beta/download>`__.