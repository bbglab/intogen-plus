

GENOME = 38
ENSEMBL ?= 92
DATASETS?=../datasets/hg${GENOME}_ensembl${ENSEMBL}_$(shell date +%Y%m%d)

# Define the variables for the recursive calls
export INTOGEN_DATASETS=$(DATASETS)
export INTOGEN_ENSEMBL=$(ENSEMBL)
export INTOGEN_GENOME=$(GENOME)

# Configure BGData
BGDATA_LOCAL=$(DATASETS)/bgdata
BGDATA_OFFLINE="FALSE"
export BGDATA_LOCAL
export BGDATA_OFFLINE

.PHONY: all bgdata clean \
	boostdm cbase combination dnds mutpanning mutrate \
	fml preprocess ptms shared smregions transvar vep

all: bgdata \
	boostdm cbase combination dnds mutpanning mutrate \
	fml preprocess ptms shared smregions transvar vep

test: bgdata shared

# Ensure DATA_DIR exists before adding files
$(DATASETS):
	mkdir $(DATASETS)

# BGData package are managed by bgdata itself
bgdata: | $(DATASETS)
	@echo bgdata datasets
	bgdata get datasets/genomereference/hg${GENOME}

clean:
	rm -rf $(DATASETS)

# the recursive calls are splitted in case one depends on another

boostdm:
	$(MAKE) -C boostDM

cbase:
	$(MAKE) -C cbase

combination:
	$(MAKE) -C combination

dnds:
	$(MAKE) -C dndscv

mutpanning:
	$(MAKE) -C mutpanning

mutrate:
	$(MAKE) -C mutrate

fml:
	$(MAKE) -C oncodrivefml

preprocess:
	@echo Building datasets for preprocessing
	$(MAKE) -C preprocess

ptms:
	$(MAKE) -C ptms

shared:
	@echo Building shared datasets
	$(MAKE) -C shared

smregions:
	$(MAKE) -C smregions

transvar:
	$(MAKE) -C transvar

vep:
	$(MAKE) -C vep