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
diagonal (\ :math:`x = y`\ ) up to an uncertainty threshold of 0.1. Genes
within the uncertainty area as well as genes with
:math:`\omega_{\text{non}} < \ 1` and :math:`\omega_{\text{mis}} < \ 1`
were deemed "ambiguous".