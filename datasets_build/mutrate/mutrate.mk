# TODO where are the files coming from

mutrate_data_srcdir = ${src_datasets}/mutrate

MUTRATE_DIR = $(DATASETS)/mutrate
$(MUTRATE_DIR): | $(DATASETS)
	mkdir $@

MUTRATE_GENOME_SIGNATURE = $(MUTRATE_DIR)/signatures.cosmic.genome.tsv
$(MUTRATE_GENOME_SIGNATURE): ${mutrate_data_srcdir}/signatures.cosmic.genome.tsv | $(MUTRATE_DIR)
	cp -f $< $@

# TODO do not name it cosmic in the output
MUTRATE_EXOME_SIGNATURE = $(MUTRATE_DIR)/signatures.cosmic.exome.tsv
$(MUTRATE_EXOME_SIGNATURE): ${mutrate_data_srcdir}/cosmic2exome.py $(MUTRATE_GENOME_SIGNATURE) $$(COUNT_CDS) $$(COUNT_WG) | $(MUTRATE_DIR)
	@echo Building mutrate exome signature
	python $< $(MUTRATE_GENOME_SIGNATURE) $(COUNT_CDS) $(COUNT_WG) $@

ALL_DATASETS += $(MUTRATE_GENOME_SIGNATURE) $(MUTRATE_EXOME_SIGNATURE)