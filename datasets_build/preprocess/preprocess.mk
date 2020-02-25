
SRC_PREPROCESS = ${DATASETS_SOURCE_FOLDER}/preprocess

DATASETS_PREPROCESS = $(DATASETS)/preprocess
$(DATASETS_PREPROCESS): | $(DATASETS)
	mkdir $@

# TODO: the file is also available in bgdata. Change the point path to it
COVERAGE_FILENAME = hg${GENOME}_100bp.coverage.regions.gz
COVERAGE = $(DATASETS_PREPROCESS)/${COVERAGE_FILENAME}
$(COVERAGE): | $(DATASETS_PREPROCESS)
	cp -f ${SRC_PREPROCESS}/${COVERAGE_FILENAME} $@

LIFTOVER = $(shell awk -v genome=hg${GENOME} '{if ($$1 == genome) {print $$2}}' ${SRC_PREPROCESS}/liftover.txt)
LIFTOVER_FILE = $(addprefix $(DATASETS_PREPROCESS)/, ${LIFTOVER})
$(DATASETS_PREPROCESS)/%ToHg${GENOME}.over.chain.gz: | $(DATASETS_PREPROCESS)
	wget -c http://hgdownload.soe.ucsc.edu/goldenPath/$*/liftOver/$*ToHg${GENOME}.over.chain.gz \
		-O $@

TARGETS_DATASETS += $(COVERAGE) $(LIFTOVER_FILE)