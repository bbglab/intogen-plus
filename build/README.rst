
Build
=====

This directory contains the required code to build
the datasets and containers (as
`Singularity <https://sylabs.io/docs/>`_
images).

Although you can find many ``*.mk`` files,
they are meant to be called form the Makefile in this directory
and not directly, as paths, or variables might not be properly defined.
However, you can get an idea of what is the build process for each dataset.


Datasets
--------

To build the datasets, you can use the Makefile in this
directory to recursively call all the Makefiles.

Currently there are 3 variables that can be defined
``ENSEMBL`` refers to the ensembl version
ENSEMBL ?= "92"
ENSEMBL_ARCHIVE ?= "apr2018"
DATASETS?=../datasets/hg$(GENOME)_ensembl$(ENSEMBL)_$(shell date +%Y%m%d)


CORES


Requirements
************

Requires a working Internet connection
and the following software:

awk
curl
make
mysql
singularity (version 2.6.1)
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
