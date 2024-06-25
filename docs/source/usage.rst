
Installation
------------

The IntOGen pipeline requires `Nextflow <https://www.nextflow.io/>`_
and `Singularity <https://sylabs.io/docs/>`_ in order to run.

Beside them, a number of different datasets need to be downloaded,
and other pieces of software installed (as Singularity containers).

For information on how to download and build all these requirements,
check the `README <https://github.com/bbglab/intogen-plus/blob/v2024/build/>`_
file in the build folder.


Usage
^^^^^

Once all the prerequisites are available, running the pipeline
only requires to execute the :file:`intogen.nf` file with the appropiate
parameters. E.g.:

.. code-block:: bash

      nextflow run intogen.nf --input <input>

There are a number of parameters and options that can be added:


-resume  Nextflow feature to allow for resumable executions.

--output <path>   Path of the input. See below for more details.

--container <path>   Path where to store the output. Default: ``intogen_analysis``.

--datasets <path>   Path to the folder containing the singularity images. Default: ``containers``.

--annotations <file>    Path to the default annotations file. Default: ``config/annotations.txt``. See the input section for more details.

--seed <int>   Seed to be used for reproducibility. This applies to 4 methods: smRegions, OncodriveCLUSTL, OncodriveFML, dNdScv.

--debug <bool>    Ask methods for a more verbose output if set to True.


Input & output
^^^^^^^^^^^^^^

Input
*****

Although the pipeline does most of its computations at the cohort level,
the pipeline is prepared to work with multiple cohorts at the same time.

Each cohort must contain, at least, the chromosome, position, ref, alt
and sample. Files are expected to be TSV files with a header line.

.. important:: All mutations should be mapped to the positive strand.
   The strand value is ignored.

In addition, each cohort must be associated with:

- cohort ID (``DATASET``): a unique identifier for each cohort.
- a cancer type (``CANCER``): although any acronym can be used here, we
  recommend to restrict to the acronyms that can be found
  in :file:`extra/data/dictionary_long_name.json`.
- a sequencing platform (``PLATFORM``): ``WXS`` for whole exome sequencing
  and ``WGS`` for whole genome sequencing
- a reference genome (``GENOMEREF``): only ``HG38`` and ``HG19`` are supported

Cohort file names, as well as the fields mentioned above
must not contain dots.

The way to provide those values is through `OpenVariant <https://github.com/bbglab/openvariant>`__ , 
a comprehensive Python package that provides different functionalities to read, parse and operate 
different multiple input file formats (e. g. tsv, csv, vcf, maf, bed). 
Whether you are planning to run single or multiple cohorts, you would 
need to provide an annotation file in yaml format to specify the above mentioned structure required by IntOGen. 
Instructions on how to build an annotation file are documented here: `OpenVariant annotation file <https://openvariant.readthedocs.io/en/latest/user_guide/annotation_structure.html>`__ .


Output
******

By default this pipeline outputs 4 files:

- :file:`cohorts.tsv`: summary of the cohorts that have been analyzed
- :file:`drivers.tsv`: summary of the results of the driver discovery by cohort
- :file:`mutations.tsv`: summary of all the mutations analyzed by cohort
- :file:`unique_drivers.tsv`: information on the genes reported as drivers (in any cohort)
- :file:`unfiltered_drivers.tsv`: information on the filters applied to the post-processing step: from the output of the combination to the final set of driver genes.

Those files can be found in the path indicated with the
``--output`` options.

Moreover, the ``--debug true`` options will generate a
:file:`debug` folder under the output folder, in which
all the input and output files of the different methods are
linked.
