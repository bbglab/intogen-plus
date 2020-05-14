# TODO where are the files coming from

mutrate_data_srcdir = ${src_datasets}/mutrate

mutrate_dir = $(INTOGEN_DATASETS)/mutrate

$(mutrate_dir): | $(INTOGEN_DATASETS)
	mkdir $@


MUTRATE_GENOME_SIGNATURE = $(mutrate_dir)/signatures.cosmic.genome.tsv
$(MUTRATE_GENOME_SIGNATURE): ${mutrate_data_srcdir}/signatures.cosmic.genome.tsv | $(mutrate_dir)
	cp -f $< $@


# TODO do not name it cosmic in the output
MUTRATE_EXOME_SIGNATURE = $(mutrate_dir)/signatures.cosmic.exome.tsv
$(MUTRATE_EXOME_SIGNATURE): ${mutrate_data_srcdir}/cosmic2exome.py $(MUTRATE_GENOME_SIGNATURE) $$(COUNT_CDS) $$(COUNT_WG) | $(mutrate_dir)
	@echo Building mutrate exome signature
	python $< $(MUTRATE_GENOME_SIGNATURE) $(COUNT_CDS) $(COUNT_WG) $@

DATASETS_TARGETS += $(MUTRATE_GENOME_SIGNATURE) $(MUTRATE_EXOME_SIGNATURE)