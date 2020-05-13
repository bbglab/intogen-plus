SHELL := /bin/bash

.PHONY: all clean
# trick to make all the first target
# and set it to the end once all variables are defined
all: containers datasets

# hg${GENOME}_ensembl${ENSEMBL}_cadd${CADD}
# TODO set the other paths
#build_dir = ../../intogen_build
#DATASETS ?= ${build_dir}/datasets
#CONTAINERS ?= ${build_dir}/containers

build_dir = ../
DATASETS ?= ${build_dir}/datasets/hg38_vep92_v0.1
CONTAINERS ?= ${build_dir}/containers/hg38_vep92_v0.1

# hg version
genome = 38
# TODO remove ensemble limitation due to BGVEP
ensembl ?= 92
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
DATASETS_TARGETS =
CONTAINERS_SUDO =
CONTAINERS_USER =

# useful variables
ensembl_db = "homo_sapiens_core_${ensembl}_${genome}"
ensembl_archive = $(shell awk -v release=${ensembl} '{if ($$1 == release) {print $$2}}' ensembl_archive.txt)
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
ENSEMBL = $(DATASETS)/.ensembl${ensembl}
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

datasets: $(DATASETS_TARGETS) | $(DATASETS)

containers: $(CONTAINERS_USER) $(CONTAINERS_SUDO) | $(CONTAINERS)

sudo: $(CONTAINERS_SUDO) | $(CONTAINERS)

clean:
	rm -rf $(DATASETS) $(CONTAINERS)
