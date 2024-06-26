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
were flagged as 'No driver' in the DRIVER column at the unfiltered_drivers.tsv file.

Additionally, candidate driver genes were further filtered using the
following criteria:

1. We discarded non-expressed genes using TCGA expression data. For tumor types directly mapping to cohorts from TCGA --including TCGA cohorts-- we removed non-expressed genes in that tumor type. We used the following criterion for non-expressed genes: genes where at least 80% of the samples showed a negative log2 RSEM. For those tumor types which could not be mapped to TCGA cohorts this filtering step was not done.
2. We also discarded genes highly tolerant to Single Nucleotide Polymorphisms (SNP) across human populations. Such genes are more susceptible to calling errors and should be taken cautiously. More specifically, we downloaded transcript specific constraints from gnomAD (release 2.1; 2018/02/14) and used the observed-to-expected ratio score (oe) of missense (mys), synonymous (syn) and loss-of-function (lof) variants to detect genes highly tolerant to SNPs. Genes enriched in SNPs (oe_mys > 1.5 or oe_lof > 1.5 or oe_syn > 1.5) with a number of mutations per sample greater than 1 were discarded. Additionally, we discarded mutations overlapping with germline variants (germline count > 5) from a panel of normals (PON) from Hartwig Medical Foundation (\ https://storage.googleapis.com/hmf-public/HMFtools-Resources/dna_pipeline/v5_31/38/variants/SageGermlinePon.98x.38.tsv.gz \ ).
3. We also discarded genes that are likely false positives according to their known function from the literature. We convened that the following genes are likely false positives: i) known long genes such as TTN, OBSCN, RYR2, etc.; ii) olfactory receptors from HORDE (\ http://bioportal.weizmann.ac.il/HORDE/\ ; Build #44c - 30 July 2019 ); iii) genes not belonging to Tier1 CGC genes lacking literature references according to CancerMine [2]_ (\ http://bionlp.bcgsc.ca/cancermine/\ ; As of 7 December 2021).
4. We also removed non CGC genes with more than 2 mutations in one sample. This abnormally high number of mutations in a sample may be the result of either a local hypermutation process or cross contamination from germline variants.
5. Finally we discarded genes whose mutations are likely the result of local hypermutation activity. More specifically, some coding regions might be the target of mutations associated to COSMIC Signature 9 (\ https://cancer.sanger.ac.uk/cosmic/signatures\) which is associated to non-canonical AID activity in lymphoid tumours. In those cancer types were Signature 9 is known to play a significant mutagenic role (i.e., LYMPHOID, CLLSLL, DLBCLNOS, NHL, ALL, LNM, TLL, BLL and HL) we discarded genes where more than 50% of mutations in a cohort of patients were associated with Signature 9.

Candidate driver genes that were not discarded composed the catalog of driver genes.

Classification according to MSKCC oncotree
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We then annotated the catalog of highly confident driver genes according
to their annotation level in CGC [1]_. Specifically, we created a three-level
annotation: i) the first level included driver genes with a reported
involvement in the source tumor type according to the CGC; ii) the
second group included CGC genes lacking reported association with the
tumor type; iii) the third group included genes that were not present in
CGC.

To match the tumor type of our analyzed cohorts and the nomenclature/acronyms of cancer types reported in the CGC we used MSKCC oncotree (as of November 2021). Resulting in 889 cancer type nodes. We customized the oncotree according to the following rules: 

1. NON_SOLID node added after TISSUE and before MYELOID and LYMPHOID
2. SOLID node added after TISSUE and before the rest of tissues
3. ALL node added after LMN and before BLL and TLL

.. note:: 
    The current version of the oncotree used in IntOGen 2023 is available at this GitHub repo: `bbglab/oncotree <https://github.com/bbglab/oncotree>`__ .

Mode of action of driver genes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We computed the mode of action for highly confident driver genes. To do
so, we first applied the aforementioned algorithm (see Mode of
action section in mutational features for more details on how the algorithm determines
the role of driver genes according to their distribution of mutations in
a cohort of samples) to classify driver genes into the three possible
roles: Act (activating or oncogene), LoF (loss-of-function or tumor
suppressor) or Amb (ambiguous or non-defined). We then combined these
predictions with prior knowledge from the Cancer Genome Census (CGC)
[3]_ using the following approach: if the inferred mode of action matched the CGC
mode of action or if the CGC mode of action was "Unknown", we used the 
inferred mode of action; otherwise, we used the CGC mode of action.

.. [1] Sondka Z, et al. The COSMIC Cancer Gene Census: describing genetic dysfunction across all human cancers. Nat Rev Cancer. 2018;18(11):696–705. doi:10.1038/s41568-018-0060-1
.. [2] Lever J, et al. CancerMine: a literature-mined resource for drivers, oncogenes and tumor suppressors in cancer. Nat Methods. 2019 Jun;16(6):505-507. doi: 10.1038/s41592-019-0422-y. Epub 2019 May 20.
.. [3] See footnote 1