SHELL := /bin/bash

.PHONY: all clean
# trick to make all the first target
# and set it to the end once all variables are defined
all: datasets

build_dir = ..
INTOGEN_DATASETS ?= ${build_dir}/datasets

# hg version
genome = 38
ensembl ?= 101
cadd ?= 1.6
# number of cores to use in steps that allow them
cores ?= 1

src_datasets = datasets

# accumulate all targets
DATASETS =

# useful variables
ensembl_db = "homo_sapiens_core_${ensembl}_${genome}"
ensembl_archive = $(shell awk -v release=${ensembl} '{if ($$1 == release) {print $$2}}' ensembl_archive.txt)
biomart_url = http://${ensembl_archive}.archive.ensembl.org/biomart/martservice
grch = GRCh${genome}

# Ensure directories exists before adding files
$(INTOGEN_DATASETS):
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

#########################
.PHONY: datasets

datasets: bgdata $(DATASETS) | $(INTOGEN_DATASETS)

clean:
	rm -rf $(INTOGEN_DATASETS)
