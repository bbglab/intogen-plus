SHELL := /bin/bash

.PHONY: all clean
# trick to make all the first target
# and set it to the end once all variables are defined
all: containers datasets

build_dir = ..
INTOGEN_DATASETS ?= ${build_dir}/datasets
INTOGEN_CONTAINERS ?= ${build_dir}/containers


# hg version
genome = 38
# TODO remove ensemble limitation due to BGVEP
ensembl ?= 92
# TODO upgrade CADD version
cadd ?= 1.4
# number of cores to use in steps that allow them
# FIXME set cores to 1
cores ?= 54


src_datasets = datasets
src_containers = containers

container_builder = ${src_containers}/build.sh


# accumulate all targets
DATASETS =
CONTAINERS_SUDO =
CONTAINERS_USER =

# useful variables
ensembl_db = "homo_sapiens_core_${ensembl}_${genome}"
ensembl_archive = $(shell awk -v release=${ensembl} '{if ($$1 == release) {print $$2}}' ensembl_archive.txt)
biomart_url = http://${ensembl_archive}.archive.ensembl.org/biomart/martservice
grch = GRCh${genome}

# Ensure directories exists before adding files
$(INTOGEN_DATASETS):
	mkdir -p $@

$(INTOGEN_CONTAINERS):
	mkdir -p $@

# Create checkpoints files so that if versions are changed
# proper files are rebuild
GENOME = $(INTOGEN_DATASETS)/.hg${genome}
ENSEMBL = $(INTOGEN_DATASETS)/.ensembl${ensembl}
CADD = $(INTOGEN_DATASETS)/.cadd${cadd}

$(GENOME): | $(INTOGEN_DATASETS)
	touch $@
$(ENSEMBL): | $(INTOGEN_DATASETS)
	touch $@
$(CADD): | $(INTOGEN_DATASETS)
	touch $@



# Use second expansion for mixed dependencies
.SECONDEXPANSION:

include ${src_datasets}/*.mk
include ${src_datasets}/*/*.mk
include ${src_containers}/*/*.mk

#########################
.PHONY: datasets containers sudo

datasets: $(DATASETS) | $(INTOGEN_DATASETS)

containers: $(CONTAINERS_USER) $(CONTAINERS_SUDO) | $(INTOGEN_CONTAINERS)

sudo: $(CONTAINERS_SUDO) | $(INTOGEN_CONTAINERS)

clean:
	rm -rf $(INTOGEN_DATASETS) $(INTOGEN_CONTAINERS)
