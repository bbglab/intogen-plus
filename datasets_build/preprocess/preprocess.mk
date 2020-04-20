
preprocess_data_srcdir = ${src_datasets}/preprocess


preprocess_dir = $(DATASETS)/preprocess
$(preprocess_dir): | $(DATASETS)
	mkdir $@


# TODO: the file is also available in bgdata. Change the point path to it
coverage_filename = hg${genome}_100bp.coverage.regions.gz
coverage_file = ${preprocess_data_srcdir}/${coverage_filename}
COVERAGE = $(preprocess_dir)/${coverage_filename}

$(COVERAGE): $(coverage_file) $$(GENOME) | $(preprocess_dir)
	cp -f $(coverage_file) $@


# TODO remove (this is just here to be able to reproduce the current version, but I think is not needed)
COVERAGE2 = $(preprocess_dir)/hg19_100bp.coverage.regions.gz

$(COVERAGE2): | $(preprocess_dir)
	cp -f ${preprocess_data_srcdir}/hg19_100bp.coverage.regions.gz $@


liftover_38_19 = hg38ToHg19.over.chain.gz
LIFTOVER38TO19 = $(preprocess_dir)/${liftover_38_19}

$(LIFTOVER38TO19): | $(preprocess_dir)
	wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz \
		-O $@


liftover_19_38 = hg19ToHg38.over.chain.gz
LIFTOVER19TO38 = $(preprocess_dir)/${liftover_19_38}

$(LIFTOVER19TO38): | $(preprocess_dir)
	wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz \
		-O $@

DATASETS_COMPUTE += $(COVERAGE) $(COVERAGE2)
DATASETS_DOWNLOAD += $(LIFTOVER19TO38) $(LIFTOVER38TO19)