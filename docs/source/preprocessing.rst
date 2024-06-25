Preprocessing
--------------

Given the heterogeneity of the datasets analyzed in the current release
of intOGen (resulting from e.g. differences in the genome aligners,
variant calling algorithms, sequencing coverage, sequencing strategy),
we implemented a pre-processing strategy aiming at reducing
possible biases. Specifically, we conducted the following filtering
steps:

1. The pipeline is configured to run using GRCh38 as reference genome. Therefore, for each input dataset the pipeline requires that the reference genome is defined. Datasets using GRCh37 as reference genome were lifted over using PyLiftover (`https://pypi.org/project/pyliftover <https://pypi.org/project/pyliftover/>`__/; version 0.4) to GRCh38. Mutations failing to liftover from GRCh37 to GRCh38 were discarded.

2. We removed mutations with equal alternate and reference alleles, duplicated mutations within the sample sample, mutations with ‘N’ as reference or alternative allele, mutations with a reference allele not matching its reference genome and mutations within non-canonical chromosomes (i.e., mutations outside chr1 to chr22, chrX and chrY).

3. Additionally, we removed mutations with low pileup mappability, i.e. mutations in regions that could potentially map elsewhere in the genome. For each position of the genome we computed the pileup mappability, defined as the average uniqueness of all the possible reads of 100bp overlapping a position and allowing up to 2 mismatches. This value is equal to 1 if all the reads overlapping a mutation are uniquely mappable while it is close to 0 if most mapping reads can map elsewhere in the genome. Positions with a pileup mappability lower than 0.9 were removed from further analyses.

4. We filtered out multiple samples from the same donor. The analysis of positive selection in tumors requires that each sample in a cohort is independent from the other samples. That implies that if the input dataset includes multiple samples from the same patient --resulting from different biopsy sites, time points or sequencing strategies-- the pipeline automatically selects the first according to its alphabetical order. Therefore, all mutations in the discarded samples are not considered anymore.

5. We also filtered out hypermutated samples. Samples carrying more than 1000 mutations for WXS and 10000 for WGS and a mutation count greater than 1.5 times the interquartile range length above the third quartile in their respective dataset were considered hypermutated and therefore removed from further analyses.

6. Datasets with filtered synonymous variants are not runnable. Most cancer driver identification methods need synonymous variants to fit a background mutation model. Therefore, datasets with less than 5 synonymous and datasets with a missense/synonymous ratio greater than 10 were excluded .

7. When the Variant Effect Predictor [1]_ (VEP) mapped one mutation into multiple transcripts associated with HUGO symbols, we selected the MANE transcript of the first HUGO symbol in alphabetical order.

8. We also discarded mutations mapping into genes without MANE transcript in VEP.111.

.. [1] McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GR, Thormann A, Flicek P, Cunningham F. The Ensembl Variant Effect Predictor. Genome Biology Jun 6;17(1):122. (2016) doi:10.1186/s13059-016-0974-4