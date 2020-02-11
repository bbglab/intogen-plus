SHELL := /bin/bash

GENOME = 38
ENSEMBL ?= 92
# TODO upgrade CADD version
CADD ?= 1.4
DATASETS?=../datasets/hg${GENOME}_ensembl${ENSEMBL}_$(shell date +%Y%m%d)
CORES ?= 1

# TODO
# export COSMIC_KEY=$(echo "email@example.com:mycosmicpassword" | base64)

.PHONY: all alltargets clean

# accumulate all targets
ALL_TARGETS=

# useful variables
ENSEMBL_DATABASE="homo_sapiens_core_${ENSEMBL}_${GENOME}"
ENSEMBL_ARCHIVE = $(shell awk -v release=${ENSEMBL} '{if ($$1 == release) {print $$2}}' ensembl_archive.txt)
BIOMART_URL = http://${ENSEMBL_ARCHIVE}.archive.ensembl.org/biomart/martservice

# trick to make all the first target
# and set it to the end once all variables are defined
all: alltargets



# The include order is important, as variables must be defined
# before being used
# The name of the variables is also important. They should not collide
include bgdata/bgdata.mk

include cbase/cbase.mk
include cgc/cgc.mk
include mutpanning/mutpanning.mk
include preprocess/preprocess.mk
include shared/shared.mk


include oncodrivefml/oncodrivefml.mk  # requies shared

# with missing parts
include smregions/smregions.mk
include transvar/transvar.mk
# dnds vep


#$(info ${ALL_TARGETS})

# Ensure DATA_DIR exists before adding files
$(DATASETS):
	mkdir $@






####### DNDSCV #######
DNDSCV_FOLDER ?= $(DATASETS)/dndscv
$(DNDSCV_FOLDER):
	mkdir $@

REF_RDA = $(DNDSCV_FOLDER)/RefCDS.rda
$(REF_RDA): $(GENOME_FASTA) | $(DNDSCV_FOLDER)
	@echo Building dNdSCV reference

# TODO this guy requires the dnds image





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

alltargets: $(ALL_TARGETS) | $(DATASETS)
#	@echo $(ALL_TARGETS)

clean:
	rm -rf $(DATASETS)
