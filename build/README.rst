
Build
=====

This directory contains a Makefile to build
the datasets and containers (as
`Singularity <https://sylabs.io/docs/>`_
images).

Although you can find many ``*.mk`` files,
they are meant to be called form the Makefile in this directory
and not directly, as paths, or variables might not be properly defined.
However, you can get an idea of what is the build process for each dataset
or container.

Currently there are several variables that can be defined:

- ``INTOGEN_DATASETS``: path where to store the datasets
- ``INTOGEN_CONTAINERS``: path where to store the containers
- ``ensembl``: specifies the ensembl version
- ``cadd``: specifies the CADD version (used for OncodriveFML)
- ``cores``: amount of cores to use by processes that can be parallelized

.. important:: Not all versions of ``ensembl`` and ``cadd``
   might work. At least, they need to be compatible with the working reference
   genome (currently hg38).

Moreover, there is a special target (``sudo``) that
can be used to build the containers that require superuser privileges
(singularity containers build by recipe).


Important notes
***************

The jar file with the latest version of MutPanning needs
to be manually downloaded from http://www.cancer-genes.org/
and placed in ``containers/mutpanning``.

The scores used by OncodriveFML are build querying directly the
CADD scores from https://cadd.gs.washington.edu/download
The process can be significantly faster and less error prone
if you download it first and replace the ``CADD_URL`` variable
in ``datasets/oncodriverfml/fml.mk`` with the full path where
you have downloaded the CADD scores.

We are providing two files for HotMAPS to avoid recomputing them;
however we suggest you to recompute them. These files are:

- ``fully_described_pdb_info.txt``: generate it with ``make annotateStructures``
  as described in the `HotMAPS wiki <https://github.com/KarchinLab/HotMAPS/wiki>`_
- ``coordinates.txt.gz``: there are 2 ``HOTMAPS_COORDINATES`` targets
  in the ``datasets/hotmaps/hotmaps.mk``. Uncomment the commented lines
  and comment the uncommented ones for the same target.


Requirements
************

Requires a working Internet connection
and the following software:

awk
cut
xz
curl
make
mysql
sqlite3
singularity
tabix
python
	bgdata
	bgsignature
	bgreference
	bgvep
	click
	numpy
	pandas
	tqdm

This software (except singularity) can be installed with
`conda <https://docs.conda.io/en/latest/>`_.

We have tested it with Singularity version 2.6.1

