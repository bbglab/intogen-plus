# TODO where are the files coming from

SRC_MUTRATE = ${DATASETS_SOURCE_FOLDER}/mutrate

DATASETS_MUTRATE = $(DATASETS)/mutrate
$(DATASETS_MUTRATE): | $(DATASETS)
	mkdir $@

MUTRATE_GENOME_SIGNATURE = $(DATASETS_MUTRATE)/signatures.cosmic.genome.tsv
$(MUTRATE_GENOME_SIGNATURE): ${SRC_MUTRATE}/signatures.cosmic.genome.tsv | $(DATASETS_MUTRATE)
	cp -f $< $@

# TODO do not name it cosmic in the output
MUTRATE_EXOME_SIGNATURE = $(DATASETS_MUTRATE)/signatures.cosmic.exome.tsv
$(MUTRATE_EXOME_SIGNATURE): $(MUTRATE_GENOME_SIGNATURE) $$(COUNT_CDS) $$(COUNT_WG) ${SRC_MUTRATE}/cosmic2exome.py | $(FOLDER)
	@echo Building mutrate exome signature
	python ${SRC_MUTRATE}/cosmic2exome.py $(MUTRATE_GENOME_SIGNATURE) $(COUNT_CDS) $(COUNT_WG) $@

TARGETS_DATASETS += $(MUTRATE_GENOME_SIGNATURE) $(MUTRATE_EXOME_SIGNATURE)