
SRC_PREPROCESS = ${DATASETS_SOURCE_FOLDER}/preprocess

DATASETS_PREPROCESS = $(DATASETS)/preprocess
$(DATASETS_PREPROCESS): | $(DATASETS)
	mkdir $@

# TODO: the file is also available in bgdata. Change the point path to it
COVERAGE_FILENAME = hg${GENOME}_100bp.coverage.regions.gz
COVERAGE = $(DATASETS_PREPROCESS)/${COVERAGE_FILENAME}
$(COVERAGE): | $(DATASETS_PREPROCESS)
	cp -f ${SRC_PREPROCESS}/${COVERAGE_FILENAME} $@

# TODO remove (this is just here to be able to reproduce the current version, but I think is not needed)
COVERAGE2 = $(DATASETS_PREPROCESS)/hg19_100bp.coverage.regions.gz
$(COVERAGE2): | $(DATASETS_PREPROCESS)
	cp -f ${SRC_PREPROCESS}/hg19_100bp.coverage.regions.gz $@

#LIFTOVER = $(shell awk -v genome=hg${GENOME} '{if ($$1 == genome) {print $$2}}' ${SRC_PREPROCESS}/liftover.txt)
#LIFTOVER_FILE = $(addprefix $(DATASETS_PREPROCESS)/, ${LIFTOVER})
#$(DATASETS_PREPROCESS)/%ToHg${GENOME}.over.chain.gz: | $(DATASETS_PREPROCESS)
#	wget -c http://hgdownload.soe.ucsc.edu/goldenPath/$*/liftOver/$*ToHg${GENOME}.over.chain.gz \
#		-O $@

LIFTOVER38TO19 = hg38ToHg19.over.chain.gz
LIFTOVER38TO19_FILE = $(addprefix $(DATASETS_PREPROCESS)/, ${LIFTOVER38TO19})
$(LIFTOVER38TO19_FILE): | $(DATASETS_PREPROCESS)
	wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz \
		-O $@

LIFTOVER19TO38 = hg19ToHg38.over.chain.gz
LIFTOVER19TO38_FILE = $(addprefix $(DATASETS_PREPROCESS)/, ${LIFTOVER19TO38})
$(LIFTOVER19TO38_FILE): | $(DATASETS_PREPROCESS)
	wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz \
		-O $@

LIFTOVER_FILE = $(LIFTOVER38TO19_FILE) $(LIFTOVER19TO38_FILE)

TARGETS_DATASETS += $(COVERAGE) $(LIFTOVER_FILE) $(COVERAGE2)