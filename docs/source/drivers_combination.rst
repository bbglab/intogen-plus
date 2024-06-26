Combining the outputs of driver identification methods
------------------------------------------------------

Rationale
^^^^^^^^^

Our goal is to provide a catalogue of driver elements which
appropriately reflects the consensus from the DIMs we run.

To combine the results of individual statistical tests, p-value
combination methods continue to be a standard approach in the field:
e.g., Fisher [1]_, Brown [2]_, [3]_ and Stouffer Z-score [4]_ methods have
been used for this purpose. These methods are useful for combining
probabilities in meta-analysis, hence to provide a ranking based on
combined significance under statistical grounds. However, the
application of these methods may bear some caveats:

1. The ranking resulting from p-value combination may lead to inconsistencies when compared to the individual rankings, i.e., they may yield a consensus ranking that does not preserve recurrent precedence relationships found in the individual rankings.

2. Some methods, like Fisher’s or Brown’s method, tend to bear anti-conservative performance, thus leading to many likely false discoveries.

3. Balanced (non-weighted) p-value combination methods may lead to faulty results just because of the influence of one or more DIM performing poorly for a given dataset.

Weighted methods to combine p-values, like the weighted Stouffer
Z-score, provide some extra room for proper balancing, in the sense of
incorporating the relative credibility of each DIM. We reasoned that any
good operational criteria to allocate weights should satisfy the
following requirements: i) provide weighting on a cohort-specific basis,
thereby allowing the relative credibility of a DIM to depend on the
cohort; ii) reflect prior knowledge about known bona-fide driver genes;
iii) reflect prior knowledge about the criteria that each DIM employed
to yield its output.

Our approach works independently for each cohort: to create a consensus
list of driver genes for each cohort, we first determine how credible
each DIM is when applied to this specific cohort, based on how
many bona-fide cancer genes reported in the COSMIC Cancer Gene Census
database (`CGC v99 <https://cancer.sanger.ac.uk/census>`__) are highly ranked according to the DIM. Once a
credibility score has been established, we use a weighted method for combining
the p-values that each DIM gives for each candidate gene: this
combination takes the DIMs credibility into account. Based on the
combined p-values, we conduct FDR correction to conclude a ranking of
candidate driver genes alongside q-values.

|image1|

Weight Estimation by Voting
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The relative credibility for each method is based on the ability of the
method to give precedence to well-known genes already collected in the
CGC catalogue of driver genes. As each method yields a ranking of driver
genes, these lists can be combined using a voting system --Schulze’s voting method.
The method allows us to consider each method as a voter with some voting rights
(weighting) which casts ballots containing a list of candidates sorted by precedence.
Schulze’s method takes information about precedence from each individual method and
produces a new consensus ranking [5]_.

Instead of conducting balanced voting, we tune the voting rights of the
methods so that we the enrichment of CGC genes at the top
positions of the consensus list is maximized. We limit the share each
method can attain in the credibility simplex --up to a uniform threshold.
The resulting voting rights are deemed the relative credibility for each method.

Ranking Score
^^^^^^^^^^^^^

Upon selection of a catalogue of bona-fide known driver elements (CGC
catalogue of driver genes) we can provide a score for each ranking
:math:`R` of genes as follows:

:math:`E(R)\  = \sum_{i=1}^N \frac{p_{i}}{\log(i + 1)}`

where :math:`p_{i}` is the proportion of elements with rank higher than
:math:`i` which belong to CGC and N is a suitable threshold to consider
only the N top ranked elements. Using :math:`E` we can define a function
:math:`f` that maps each weighting vector :math:`w` (in the 5-simplex)
to a value :math:`E(R_{w})` where :math:`R_{w}` denotes the consensus
ranking obtained by applying Schulze’s voting with voting rights given
by the weighting vector :math:`w`.

Optimization with constraints
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally we are bound to find a good candidate for
:math:`\widehat{w} = \textrm{argmax}(f)`. For each method to have chances
to contribute in the consensus score, we impose the mild constraint of
limiting the share of each method up to 0.3.

Optimization is then carried out in two steps: we first find a good
candidate :math:`\widehat{w_{0}}` by exhaustive search in a rectangular grid
satisfying the constraints defined above (with grid step=0.05); in the second step
we take :math:`\widehat{w_{0}}` as the seed for a stochastic hill-climbing
procedure (we resort to Python’s scipy.optimize “basinhopping”, method=SLSQP
and stepsize=0.05).

Estimation of combined p-values using weighted Stouffer Z-score
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Using the relative weight estimate that yields a maximum value of the
objective function f we combined the p-values resorting to the weighted
Stouffer Z-score method. Thereafter we performed Benjamini-Hochberg FDR
correction with the resulting combined p-values, yielding one q-value
for each genomic element. If the element belongs to CGC, we computed its
q-value using only the collection of p-values computed for CGC genes.
Otherwise, we computed the q-value using all the computed p-values.


Tiers of driver genes from sorted list of combined rankings and p-values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To finalize the analysis we considered only genes with at least two
mutated samples in the cohort under analysis. These genes were
classified into four groups according to the level of evidence in that
cohort that the gene harbours positive selection.

1) The first group, named as TIER1, contained genes showing high confidence and agreement in their positive selection signals. Given the ranked list of genes obtained by the Schulze voting, TIER1 comprises all the ranked genes whose ranking is higher than the first gene with combined q-value lower than a specific threshold (by default threshold=0.05). 
2) The second group, name as TIER2, was devised to contain known cancer driver genes, showing mild signals of positive selection, that were not included in TIER1. More in detail, we defined TIER2 genes as those CGC genes, not included in TIER1, whose CGC q-value was lower than a given threshold (default CGC q-value=0.25). CGC q-value is computed by performing multiple test correction of combined p-values restricted to CGC genes. 
3) The third group, are genes not included in TIER1 or TIER2 with scattered signals of positive selection, frequently coming from one single method. Particularly, given the ranked list of genes by the Schulze voting, TIER3 was composed of all the ranked genes with q-value lower than a given threshold (by default threshold=0.05) whose ranking is higher than TIER1 last gene position and lower than the rejection ranking position. The rejection ranking position is defined as the ranking position for which all elements have a q-value lower than the input threshold (by default threshold=0.05). Finally, other genes not included in the aforementioned classes are considered non-driver genes.

Combination benchmark
^^^^^^^^^^^^^^^^^^^^^
.. important:: 
   This benchmark was performed on `IntOGen Plus v2020 <https://intogen.readthedocs.io/en/v2020>`__ .


We have assessed the performance of the combination compared to i) each
of the six individual methods and ii) other strategies to combine the
output from cancer driver identification methods.

Datasets for evaluation
~~~~~~~~~~~~~~~~~~~~~~~

To ensure a reliable evaluation we decided to perform an evaluation
based on the 32 Whole-Exome cohorts of the TCGA PanCanAtlas project
(downloaded from
`*https://gdc.cancer.gov/about-data/publications/pancanatlas* <https://gdc.cancer.gov/about-data/publications/pancanatlas>`__).
These cohorts sequence coverage, sequence alignment and somatic mutation
calling were performed using the same methodology guaranteeing that
biases due to technological and methodological artifacts are minimal.

The Cancer Genes Census --version v87-- was downloaded from the COSMIC
data portal
(`*https://cancer.sanger.ac.uk/census* <https://cancer.sanger.ac.uk/census>`__)
and used as a positive set of known cancer driver genes.

We created a catalog of genes that are known not to be involved in
cancerogenesis. This set includes very long genes (HMCN1, TTN, OBSCN,
GPR98, RYR2 and RYR3), and a list of olfactory receptors from Human
Olfactory Receptor Data Exploratorium (HORDE)
(https://genome.weizmann.ac.il/horde/; download date 14/02/2018).
In addition, for all TCGA cohorts, we added non-expressed genes, defined
as genes where at least 80% of the samples showed a RSEM expressed in
log2 scale smaller or equal to 0. Expression data for TCGA was
downloaded from
`*https://gdc.cancer.gov/about-data/publications/pancanatlas* <https://gdc.cancer.gov/about-data/publications/pancanatlas>`__.

Metrics for evaluation
~~~~~~~~~~~~~~~~~~~~~~

We defined a metric, referred to as CGC-Score, that is intended to
measure the quality of a ranking of genes as the enrichment of CGC
elements in the top positions of the ranking; specifically given a
ranking :math:`R` mapping each element to a rank, we define the
CGC-Score of :math:`R` as

:math:`\text{CGC-Score}(R)\  = \sum_{i=1}^N\frac{p_{i}}{log(i + 1)} \; /\; \sum_{i=1}^N\frac{1}{log(i + 1)}`

where :math:`p_{i}` is the proportion of elements with rank
:math:`\leq i` that belong to CGC and :math:`N` is a suitable threshold
to consider just the top elements in the ranking (by default N=40).

We estimated the CGC-Score across TCGA cohorts for the individual
methods ranking and the combined Schulze ranking.

Similarly, we defined a metric, referred to as Negative-Score, that aims
to measure the proportion non-cancer genes among the top positions in
the ranking. Particularly, given a ranking :math:`R` mapping each
element to a rank, we define the Negative-Score of :math:`R` as:

:math:`\text{Negative-Score}(R)\  = \sum_{i=1}^N \frac{p_{i}}{log(i + 1)}\; /\; \sum_{i=1}^N \frac{1}{log(i + 1)}`

where :math:`p_{i}` is the proportion of elements with rank
:math:`\leq i` that belong to the negative set and :math:`N` is a
suitable threshold to consider just the top elements in the ranking (by
default N = 40). We estimated the Negative-Score across TCGA cohorts for
the individual methods ranking and the combined Schulze ranking.

Comparison with individual methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We compared the CGC-Score and Negative-Score of our combinatorial
selection strategy with the individual output from the six driver
discovery methods integrated in the pipeline.

As a result we observed a consistent increase in CGC-Score of the
combinatorial strategy compared to individual methods across TCGA
cohorts (see Figure below panel A-B). Similarly, we observed a consistent decrease in
Negative-Score across TCGA cohorts (see Figure below panel C). In summary, the
evaluation shows that the combinatorial strategy increases the True
Positive Rate (measured using the CGC-Score) preserving lower False
Positive Rate (measured using the Negative-Score) than the six
individual methods included in the pipeline.

Leave-one-out combination
~~~~~~~~~~~~~~~~~~~~~~~~~

We aimed to estimate the contribution of each method’s ranking to the
final ranking after Schulze's weighted combination. To tackle this
question, we measured the effect of removing one method from the
combination by, first, calculating the CGC-Score of the combination
excluding the aforementioned method and, next, obtaining its ratio with
the original combination (i.e., including all methods). This was
iteratively calculated for all method across TCGA cohorts. Methods that
positively contributed to the combined ranking quality show a ratio
below one, while methods that negatively contributed to the combined
ranking show a ratio greater than one.

We observed that across TCGA cohorts most of the methods contributed
positively (i.e., ratio above one) to the weighted combination
performance. Moreover, there is substantial variability across TCGA
cohorts in the contribution of each method to the combination
performance. In summary, all methods contributed positively to the
combinatorial performance across TCGA supporting our methodological
choice for the individual driver discovery methods (see Figure below panel E).

Comparison with other combinatorial selection methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We compared the CGC-Score and Negative-Score of our combinatorial
selection strategy against other methods frequently used employed to
produce ranking combinations, either based on ranking information --such
as Borda Count [6]_ -- or based on statistical information --such as
Fisher [1]_ or Brown [2]_, [3]_ methods. Hereto, we briefly describe
the rationale of the four methods we used to benchmark our ranking. For
the sake of compact notation, let’s denote the rank and p-value of gene
:math:`g` produced by method :math:`m_{i}` as :math:`r_{i,g}` and
:math:`p_{i,g}`, respectively.

*Borda Count:* For each ranked item :math:`g` and method :math:`m_{i},`
it assigns a score :math:`s_{i,g} = N - l_{i,g},` where :math:`N`
stands for the total number of items to rank and :math:`l_{i,g}` is
the number of items ranked below :math:`g` according to method
:math:`m_{i}`. For each item :math:`g` an overall score
:math:`s_{g}= s_{1,g} + \ldots + s_{k,g}` can then be
computed for each :math:`g,` whence a ranking is obtained by descending
sort.

*Fisher:* It is based on the p-values :math:`p_{i,g}`. For each item
:math:`g` the method produces a new combined p-value by computing the
statistic:

:math:`F_{g} = - 2\log\ p_{i, g} \sim \chi_{2k}^{2}`.

Under the null hypothesis, :math:`F_{g}` are distributed as a chi-square
with :math:`2k` degrees of freedom, whence a p-value, which in turn
yields a raking by ascending sort. Its applicability is limited by the
assumption that the methods provide independent significance tests.

*Brown:* This method overcomes the independence requirement of Fisher’s
method by modeling the dependencies between the statistical tests
produced by each method, specifically by estimating the covariance
:math:`\Omega_{i,j} = \textrm{cov}( - 2\log p_{i,g}, - 2\log p_{j,g}).`
Brown’s method [2]_ and its most recent adaptation [3]_ have been
proposed as less biased alternatives to Fisher.

We then computed the CGC-Score and Negative-Score based on the consensus
ranking from the aforementioned combinatorial methods and compared them
to our Schulze’s weighted combination ranking across all TCGA cohorts.
Our combinatorial approach met a larger enrichment in known cancer genes
for 29/32 (90%) TCGA cohorts (see Figure below panel D).

|image2|


.. rubric:: Footnotes


.. [1] Fisher R.A. (1948) figure to question 14 on combining independent tests of significance. Am. Statistician , 2, 30–31.

.. [2] Brown, M. B. 400: A Method for Combining Non-Independent, One-Sided Tests of Significance. Biometrics 31, 987 (1975). DOI: 10.2307/2529826

.. [3] William Poole, et al. Combining dependent P-values with an empirical adaptation of Brown’s method, Bioinformatics, Volume 32, Issue 17, 1 September 2016, Pages i430–i436, https://doi.org/10.1093/bioinformatics/btw438

.. [4] Zaykin, D. V. Optimally weighted Z-test is a powerful method for combining probabilities in meta-analysis. Journal of Evolutionary Biology 24, 1836–1841 (2011). doi: 10.1111/j.1420-101.2011.02297.x

.. [5] https://arxiv.org/pdf/1804.02973.pdf

.. [6] Emerson P. The original Borda count and partial voting. Soc Choice Welf (2013) 40:353–358. doi 10.1007/s00355-011-0603-9


.. |image1| image:: /_static/schema_intogen_methods.png
   :width: 7.2in
   :height: 2.3in
   :align: middle
   :scale: 100%

.. |image2| image:: /_static/benchmark.png
   :width: 9.00000in
   :height: 6in
   :align: middle




