# TODO where are the files coming from

FOLDER_MUTRATE = $(DATASETS)/mutrate
$(FOLDER_MUTRATE): | $(DATASETS)
	mkdir $@

MUTRATE_GENOME_SIGNATURE = $(FOLDER_MUTRATE)/signatures.cosmic.genome.tsv
$(MUTRATE_GENOME_SIGNATURE): mutrate/signatures.cosmic.genome.tsv | $(FOLDER_MUTRATE)
	cp -f $< $@

# TODO do not name it cosmic in the output
MUTRATE_EXOME_SIGNATURE = $(FOLDER_MUTRATE)/signatures.cosmic.exome.tsv
$(MUTRATE_EXOME_SIGNATURE): $(MUTRATE_GENOME_SIGNATURE) $$(COUNT_CDS) $$(COUNT_WG) mutrate/cosmic2exome.py | $(FOLDER)
	@echo Building mutrate exome signature
	python mutrate/cosmic2exome.py $(MUTRATE_GENOME_SIGNATURE) $(COUNT_CDS) $(COUNT_WG) $@

ALL_TARGETS += $(MUTRATE_GENOME_SIGNATURE) $(MUTRATE_EXOME_SIGNATURE)