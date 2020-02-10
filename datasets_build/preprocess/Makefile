
FOLDER ?= ${INTOGEN_DATASETS}/preprocess
GENOME ?= ${INTOGEN_GENOME}
CORES ?= 1

.PHONY: all clean

COVERAGE = $(FOLDER)/hg${GENOME}_100bp.coverage.regions.gz
LIFTOVER = $(shell awk -v genome=hg${GENOME} '{if ($$1 == genome) {print $$2}}' liftover.txt)
LIFTOVER_FILE = $(addprefix $(FOLDER)/, ${LIFTOVER})

all: $(COVERAGE) $(LIFTOVER_FILE) | check

$(FOLDER):
	mkdir $@

# TODO: the file is also available in bgdata. Change the point path to it
$(COVERAGE): | $(FOLDER)
	cp -f hg${GENOME}_100bp.coverage.regions.gz $@


$(FOLDER)/%ToHg${GENOME}.over.chain.gz: | $(FOLDER)
	wget -c http://hgdownload.soe.ucsc.edu/goldenPath/$*/liftOver/$*ToHg${GENOME}.over.chain.gz \
		-O $@


check:
ifeq ($(FOLDER), /preprocess)
	$(error FOLDER not set)
endif
ifeq (${GENOME}, )
	$(error GENOME not set)
endif