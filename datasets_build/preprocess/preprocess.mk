
FOLDER_PREPROCESS = $(DATASETS)/preprocess
$(FOLDER_PREPROCESS): | $(DATASETS)
	mkdir $@

# TODO: the file is also available in bgdata. Change the point path to it
COVERAGE_FILENAME = hg${GENOME}_100bp.coverage.regions.gz
COVERAGE = $(FOLDER_PREPROCESS)/${COVERAGE_FILENAME}
$(COVERAGE): | $(FOLDER_PREPROCESS)
	cp -f preprocess/${COVERAGE_FILENAME} $@

LIFTOVER = $(shell awk -v genome=hg${GENOME} '{if ($$1 == genome) {print $$2}}' preprocess/liftover.txt)
LIFTOVER_FILE = $(addprefix $(FOLDER_PREPROCESS)/, ${LIFTOVER})
$(FOLDER_PREPROCESS)/%ToHg${GENOME}.over.chain.gz: | $(FOLDER_PREPROCESS)
	wget -c http://hgdownload.soe.ucsc.edu/goldenPath/$*/liftOver/$*ToHg${GENOME}.over.chain.gz \
		-O $@

ALL_TARGETS += $(COVERAGE) $(LIFTOVER_FILE)