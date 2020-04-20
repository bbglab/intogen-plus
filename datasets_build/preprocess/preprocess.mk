
preprocess_data_srcdir = ${src_datasets}/preprocess

PREPROCESS_DIR = $(DATASETS)/preprocess
$(PREPROCESS_DIR): | $(DATASETS)
	mkdir $@

# TODO: the file is also available in bgdata. Change the point path to it
coverage_filename = hg${genome}_100bp.coverage.regions.gz
coverage_file = ${preprocess_data_srcdir}/${coverage_filename}
COVERAGE = $(PREPROCESS_DIR)/${coverage_filename}
$(COVERAGE): $(coverage_file) $$(GENOME) | $(PREPROCESS_DIR)
	cp -f $(coverage_file) $@

# TODO remove (this is just here to be able to reproduce the current version, but I think is not needed)
COVERAGE2 = $(PREPROCESS_DIR)/hg19_100bp.coverage.regions.gz
$(COVERAGE2): | $(PREPROCESS_DIR)
	cp -f ${preprocess_data_srcdir}/hg19_100bp.coverage.regions.gz $@

liftover_38_19 = hg38ToHg19.over.chain.gz
LIFTOVER38TO19 = $(PREPROCESS_DIR)/${liftover_38_19}
$(LIFTOVER38TO19): | $(PREPROCESS_DIR)
	wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz \
		-O $@

liftover_19_38 = hg19ToHg38.over.chain.gz
LIFTOVER19TO38 = $(PREPROCESS_DIR)/${liftover_19_38}
$(LIFTOVER19TO38): | $(PREPROCESS_DIR)
	wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz \
		-O $@


ALL_DATASETS += $(COVERAGE) $(LIFTOVER19TO38) $(LIFTOVER38TO19) $(COVERAGE2)