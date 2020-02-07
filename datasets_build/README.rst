
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

awk
curl
make
mysql
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