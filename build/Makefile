SHELL := /bin/bash

GENOME = 38
ENSEMBL ?= 92
# TODO upgrade CADD version
CADD ?= 1.4
DATASETS ?= datasets/hg${GENOME}_ensembl${ENSEMBL}_cadd${CADD}
CONTAINERS ?= containers/hg${GENOME}_ensembl${ENSEMBL}_cadd${CADD}
CORES ?= 1

# TODO use bash scripts when needed

# TODO combine datasets and containers in the folder.
# one folder hg${GENOME}_ensembl${ENSEMBL}_cadd${CADD} and inside container and datasets

# TODO create a file .hg38 that forces rebuild of all the thing that need it if we change the version
# shame thing with cadd and ensembl version

DATASETS_SOURCE_FOLDER = ../datasets_build
CONTAINERS_SOURCE_FOLDER = ../containers_build

.PHONY: all alltargets clean datasets containers sudo-containers

# accumulate all targets
TARGETS_DATASETS=
TARGETS_CONTAINERS=
TARGETS_CONTAINERS_SUDO=

# useful variables
ENSEMBL_DATABASE="homo_sapiens_core_${ENSEMBL}_${GENOME}"
ENSEMBL_ARCHIVE = $(shell awk -v release=${ENSEMBL} '{if ($$1 == release) {print $$2}}' ensembl_archive.txt)
BIOMART_URL = http://${ENSEMBL_ARCHIVE}.archive.ensembl.org/biomart/martservice

# trick to make all the first target
# and set it to the end once all variables are defined
all: alltargets


# Use second expansion for mixed dependencies
.SECONDEXPANSION:

include ${DATASETS_SOURCE_FOLDER}/*/*.mk
include ${CONTAINERS_SOURCE_FOLDER}/*/*.mk


# Ensure directories exists before adding files
$(DATASETS):
	mkdir -p $@

$(CONTAINERS):
	mkdir -p $@




####### VEP #######
VEP_FOLDER = $(DATASETS)/vep
$(VEP_FOLDER): | $(DATASETS)
	mkdir $@

# TODO what is vep_install??


### Boostdm

#include mutrate/mutrate.mk  # requires shared

####### PTMS #######
# TODO not required?

####### BOOSTDM #######
# TODO not required?






#########################


datasets: $(TARGETS_DATASETS) | $(DATASETS)

containers: $(TARGETS_CONTAINERS) $(TARGETS_CONTAINERS_SUDO) | $(CONTAINERS)

sudo: $(TARGETS_CONTAINERS_SUDO) | $(CONTAINERS)

alltargets: containers datasets
#	@echo $(ALL_TARGETS)

clean:
	rm -rf $(DATASETS) $(CONTAINERS)
