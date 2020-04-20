SHELL := /bin/bash

.PHONY: all clean
# trick to make all the first target
# and set it to the end once all variables are defined
all: containers datasets

# hg${GENOME}_ensembl${ENSEMBL}_cadd${CADD}
build_dir = ../../intogen_build
DATASETS ?= ${build_dir}/datasets
CONTAINERS ?= ${build_dir}/containers

# hg version
genome = 38
ensmebl ?= 92
# TODO upgrade CADD version
cadd ?= 1.4
# number of cores to use in steps that allow them
cores ?= 1



# TODO use bash scripts when needed

src_datasets = ../datasets_build
src_containers = ../containers_build

# tmp directory used during build
# must be an absolute path
tmpdir=`mktemp -d`

# accumulate all targets
ALL_DATASETS =
CONTAINERS_SUDO =
CONTAINERS_USER =

# useful variables
ensembl_db = "homo_sapiens_core_${ensmebl}_${genome}"
ensembl_archive = $(shell awk -v release=${ensmebl} '{if ($$1 == release) {print $$2}}' ensembl_archive.txt)
biomart_url = http://${ensembl_archive}.archive.ensembl.org/biomart/martservice
grch = GRCh${genome}

# Ensure directories exists before adding files
$(DATASETS):
	mkdir -p $@

$(CONTAINERS):
	mkdir -p $@

# Create checkpoints files so that if versions are changed
# proper files are rebuild
GENOME = $(DATASETS)/.hg${genome}
ENSEMBL = $(DATASETS)/.ensembl${ensmebl}
CADD = $(DATASETS)/.cadd${cadd}

$(GENOME): | $(DATASETS)
	touch $@
$(ENSEMBL): | $(DATASETS)
	touch $@
$(CADD): | $(DATASETS)
	touch $@



# Use second expansion for mixed dependencies
.SECONDEXPANSION:

include ${src_datasets}/*/*.mk
include ${src_containers}/*/*.mk

#########################
.PHONY: datasets containers sudo

datasets: $(ALL_DATASETS) | $(DATASETS)

containers: $(CONTAINERS_USER) $(CONTAINERS_SUDO) | $(CONTAINERS)

sudo: $(CONTAINERS_SUDO) | $(CONTAINERS)

clean:
	rm -rf $(DATASETS) $(CONTAINERS)
