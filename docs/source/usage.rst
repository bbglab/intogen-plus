
Installation
------------

The IntOGen pipeline requires `Nextflow <https://www.nextflow.io/>`_
and `Singularity <https://sylabs.io/docs/>`_ in order to run.

Beside them, a number of different datasets need to be downloaded,
and other pieces of software installed (as Singularity containers).

For information on how to download and build all these requirements,
check the `README <https://bitbucket.org/intogen/intogen-plus/src/master/build/>`_
file in the build folder.


Usage
-----

Once all the prerequisites are available, running the pipeline
only requires to execute the :file:`intogen.nf` file with the appropiate
parameters. E.g.:

.. code-block:: bash

   nextflow run intogen.nf --input <input>

There are a number of parameters and options that can be added.

.. glossary::

   **-resume**
     Nextflow feature to allow for resumable executions.

   **--input <path>**
     Path of the input. See below for more details.

   **--output <path>**
     Path where to store the output. Default: ``intogen_analysis``.

   **--containers <path>**
     Path to the folder containing the singularity images.
     Default: ``containers``.

   **--datasets <path>**
     Path to the folder the datasets from
     Default: ``containers``.

   **--annotations <file>**
     Path to the default annotations file.
     Default: ``config/annotations.txt``.
     See the input section for more details.

   **--debug true**
     Show intermediate results and do not halt the execution on error.
     This is useful when dealing with many cohorts in which some
     are likely to fail.

Input & output
--------------

Input
*****

Although the pipeline does most of its computations at the cohort level,
the pipeline is prepared to work with multiple cohorts at the same time.

Each cohort must contain, at least, the chromosome, position, ref, alt
and sample. Files are expected to be TSV files with a header line.

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

The way to provide those values is through the use of
annotations and datasource files supported by the
`bgparsers package <https://bgparsers.readthedocs.io/en/latest/usage.html#command-groupby>`_.

In the case that you are interesting in analysing a single cohort,
the cohort filename until the first dot will be used as cohort ID.
The rest of the fields can be easily provided with an annotations file.
See :file:`config/annotations.txt` for an example.

When you plan to run multiple cohorts, it is recommended to
combine the simple annotations that can be achieved with
an annotations file with the power of :file:`bginfo` files
(with ``.bginfo`` extension).
One such file (which can include others), allows to create
smart annotations and input certain files from your directory structure.
The pipeline can receive the :file:`bginfo` file or the folder where it is.
In this case, the **annotations** file (passed using ``--annotations``)
allows to make some general definitions of variables that
are used in all :file:`bginfo` files.

See examples of :file:`bginfo` in the
`bgparsers documentation <https://bgparsers.readthedocs.io/en/latest/datasource.html>`_.


Output
******

By default this pipeline outputs 4 files:

- :file:`cohorts.tsv`: summary of the cohorts that have been analyzed
- :file:`drivers.tsv`: summary of the results of the driver discovery by cohort
- :file:`mutations.tsv`: summary of all the mutations analyzed by cohort
- :file:`unique_drivers.tsv`: information on the genes reported as drivers (in any cohort)

Those files can be found in the path indicated with the
``--output`` options.

Moreover, the ``--debug true`` options will generate a
:file:`debug` folder under the output folder, in which
all the input and output files of the different methods are
linked.
