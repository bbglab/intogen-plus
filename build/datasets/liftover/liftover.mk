
# LiftOver chain files for pyliftover

liftover_data_srcdir = ${src_datasets}/liftover

liftover_dir = $(INTOGEN_DATASETS)/liftover
$(liftover_dir): | $(INTOGEN_DATASETS)
	mkdir $@

$(liftover_dir)/%.over.chain.gz: ${liftover_data_srcdir}/download.sh | $(liftover_dir)
	$< $@

DATASETS += $(liftover_dir)/hg38ToHg19.over.chain.gz \
	$(liftover_dir)/hg19ToHg38.over.chain.gz