Drivers postprocessing
----------------------

The intOGen pipeline outputs a ranked list of driver genes per input
cohort. We aimed to create a comprehensive catalog of driver genes per
tumor type from all the cohorts included in this version.

Then, we performed a filtering on automatically generated driver gene
lists per cohort. This filtering is intended to lessen artifacts from
the cohort-specific driver lists (e.g., due to errors in calling
algorithms, errors introduced by positive selection methods, local
hypermutations effects, undocumented filtering of mutations, etc.).

Therefore, we first created a collection of candidate driver genes by
selecting significant non CGC genes (i.e., q-value <0.05) with, at least
two significant bidders (i.e., two methods render them as significant);
or CGC significant genes (either q-value < 0.05 or CGC q-value < 0.25)
from individual cohorts. All genes that did not fulfill these
requirements were not longer considered.

Additionally, candidate driver genes were further filtered using the
following criteria:

1. We discarded non-expressed genes using TCGA expression data. For
       tumor types directly mapping to cohorts from TCGA --or those from
       TCGA themselves-- we removed non-expressed genes in that tumor
       type. We used the aforementioned definition of non-expressed
       genes (i.e., genes where at least 80% of the samples showed a
       RSEM expressed in log2 scale less or equal to 0). Tumor types
       lacking of mapping to TCGA cohorts did not incorporate this
       filtering step.

2. We also discarded genes highly tolerant to Single Nucleotide
       Polymorphisms (SNP) across human populations. Such genes are more
       susceptible to calling errors and should be taken cautiously.
       More specifically, we downloaded transcript specific constraints
       from gnomAD (release 2.1; 14/02/2018) and used the observed /
       expected score, also called as oe score, of missense, synonymous
       and lof variants to detect genes highly tolerant to SNPs. Genes
       enriched in SNPs (oe\_mys > 1.5 or oe\_lof > 1.5 or oe\_syn >1.5)
       with a number of mutations per sample greater than 1 were
       discarded. Additionally, we discarded mutations overlapping with
       germline variants (germline count >5) from a panel of normals
       (PON) from Hartwig Medical Foundation
       (`*https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd/download?path=%2FHMF-Pipeline-Resources&files=SOMATIC\_PON.vcf.gz* <https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd/download?path=%2FHMF-Pipeline-Resources&files=SOMATIC_PON.vcf.gz>`__).

3. We also discarded genes that are likely false positives according to
       their known function from literature. We defined as likely false
       positives to i) Known very long genes such as TTN, OBSCN, RYR2,
       etc. ii) Olfactory receptors from HORDE
       (http://bioportal.weizmann.ac.il/HORDE/; download date
       14/02/2018) iii) Non Tier1 Cancer Gene Census genes lacking of
       literature references according to CancerMine
       (`*http://bionlp.bcgsc.ca/cancermine/* <http://bionlp.bcgsc.ca/cancermine/>`__).

4. We also removed non CGC genes with more than 3 mutations in one
       sample. This abnormally high number of mutations in a sample may
       be the result of either a local hypermutation process or cross
       contamination from germline variants.

5. Finally we discarded genes whose mutations are the result of local
       hypermutation events. More specifically, samples with
       immunoglobulin gene hypermutation, harbour high number of
       mutations attributed to COSMIC Signature 9
       (`*https://cancer.sanger.ac.uk/cosmic/signatures* <https://cancer.sanger.ac.uk/cosmic/signatures>`__).
       Some coding regions might be the target Signature9 associated
       mutations, such as the cases of immunoglobulins (IGLLs). In those
       cancer types were Signature9 is considered active (i.e., AML,
       Non-Hodgkin Lymphomas, B-cell Lymphomas, CLL and Myelodysplastic
       syndromes), we discarded genes where more than 50% of mutations
       in a cohort of patients are associated to Signature9.

Candidate driver genes that were not discarded composed the catalog of
driver genes.

Classification according to their level of annotation from CGC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We then annotated the catalog of highly confident driver genes according
to their level of annotation in the CGC. Concretely, we created a
three-level annotation. The first level included driver genes with a
reported involvement in the source tumor type according to the CGC. The
second group included CGC genes lacking of reported association with the
tumor type, using CGC as a source of information. The third groups
involved genes that were not present in CGC.

To perform the association between the tumor type of our analyzed
cohorts and the nomenclature/acronyms of cancer types reported in the
CGC we manually created a mapping of all the names of tumor types from
the CGC into cancer types defined in our study following the rules:

-  All the equivalent terms for a cancer type reported in the CGC using
       the Somatic Tumor Type field (e.g. breast, breast carcinoma,
       breast cancer, etc.), were mapped into the same tumor type.

-  CGC terms with an unequivocal mapping into our cancer types were
       automatically linked (e.g., breast with BRCA).

-  CGC terms representing fine tuning classification of a more prevalent
       cancer type, which did not represent an independent cohort in our
       study; were mapped into their closest parent tumor type in our
       study (e.g., malignant melanoma of soft parts into cutaneous
       melanoma or alveolar soft part sarcoma into sarcoma).

-  Adenomas were mapped to carcinomas of the same cell type (e.g.,
       hepatic adenoma into hepatic adenocarcinoma, salivary gland
       adenoma into salivary gland adenocarcinoma, etc.).

-  CGC parent terms mapping into several tumor types from our study were
       mapped into each of the potential child tumor types. For
       instance, the term non small cell lung cancer -- or NSCLC for
       short-- was mapped into LUAD (lung adenocarcinoma) and LUSC (lung
       squamous cell carcinoma).

-  Finally, CGC terms associated with benign lesions, with unspecified
       tumor types (e.g., other, other tumor types, other CNS, etc.) or
       with tumor types with unavailable parent in our study were left
       unmatched.

Mode of action of driver genes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We computed the mode of action for highly confident driver genes. To do
so, we first performed a pan-cancer run of dNdScv across all TCGA
cohorts. We then applied the aforementioned algorithm (see Mode of
action section for more information about how the algorithm determines
the role of driver genes according to their distribution of mutations in
a cohort of samples) to classify driver genes into the three possible
roles: Act, LoF or Amb. We then combined these predictions with prior
knowledge from the Cancer Genome Interpreter following the next rules:
i) when the mode of action from our classification and prior knowledge
agreed we used the consensus mode of action; ii) when the gene was not
included in the prior knowledge list we selected the inferred mode of
action iii) when there was no agreement between methods, we selected the
mode of action from the prior knowledge list.