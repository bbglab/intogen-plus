
Build
=====

This directory contains a Makefile to build
the datasets and containers (as
`Singularity <https://sylabs.io/docs/>`_
images). To run it:

.. code:: bash

   cd build
   make

The following environment variable *+must** be defined:

- ``COSMIC_KEY``: it is used to access and retrieve the CGC
  dataset from `COSMIC <https://cancer.sanger.ac.uk/census>`_.
  It can be defined as: ``export COSMIC_KEY=$(echo "<email>:<password>" | base64)``

Currently there are several environment variables that can be defined:

- ``INTOGEN_DATASETS``: path where to store the datasets.
  Default ``../datasets``.
- ``INTOGEN_CONTAINERS``: path where to store the containers.
  Default ``../containers``.
- ``ensembl``: specifies the ensembl version
- ``cadd``: specifies the CADD version (used for OncodriveFML)
- ``cores``: amount of cores to use by processes that can be parallelized

.. important:: Not all versions of ``ensembl`` and ``cadd``
   might work. At least, they need to be compatible with the working reference
   genome (currently hg38).

Moreover, there is a special make target (``sudo``) that
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

Less important notes
********************

Although you can find many ``*.mk`` files,
they are meant to be called form the Makefile in this directory
and not directly, as paths, or variables might not be properly defined.
However, you can get an idea of what is the build process for each dataset
or container.

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
and the following software::

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
		bgparsers
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


Updating versions
*****************

By default, once it has been downloaded, any dataset that is
not specified in the ``bgdata.mk`` will not be updated.
In order to update a particular dataset, it must be removed
from the filesystem and make re-executed. Only the missing datasets
and any other dataset that depends on it will be updated.

Regarding the containers, some of them point to latest versions
of the software they install. In those cases, the container must
be manually deleted in other to be updated.

Unfortunately, not all recent version of ensembl can be used out of the box.
To be able to specify a version above 101 in the ``ensembl`` parameter of the
makefile, you need to update the ``ensembl_archive.txt`` file, which maps
the ensemble version with the release date parameter used in queries.
