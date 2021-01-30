.PHONY: all clean cohorts check-env

# The following environment variables must be defined
INTOGEN_DATASETS ?= ../datasets
export INTOGEN_DATASETS
# FOLDERS with input datasets
INPUTS = ${INTOGEN_INPUTS}
PIPELINE = ${INTOGEN_RESULTS}


CGC = ${INTOGEN_DATASETS}/cgc/cancer_gene_census_parsed.tsv
BIOMART = ${INTOGEN_DATASETS}/regions/cds_biomart.tsv
VEP = ${INTOGEN_DATASETS}/vep/vep.tsv.gz

# OUTPUTS
OUTPUT_FOLDER = output
CANCERS = $(OUTPUT_FOLDER)/ctypes.json
STATS_COHORTS = $(OUTPUT_FOLDER)/cohorts.tsv
GENES = $(OUTPUT_FOLDER)/genes.tsv
DRIVER_GENES = $(OUTPUT_FOLDER)/drivers.tsv
MUTATIONS = $(OUTPUT_FOLDER)/mutations.tsv.gz


# Other intermediate files
OUTPUT_FOLDER_DRIVERS = $(OUTPUT_FOLDER)/drivers

all: cohorts cancers genes drivers mutations

cancers: $(CANCERS)
cohorts: $(STATS_COHORTS)
genes: $(GENES)
drivers: check-env $(DRIVER_GENES)
mutations: $(MUTATIONS)


# Ensure directories exits
$(OUTPUT_FOLDER):
	mkdir $(OUTPUT_FOLDER)

$(OUTPUT_FOLDER_DRIVERS): | $(OUTPUT_FOLDER)
	mkdir $(OUTPUT_FOLDER_DRIVERS)


$(CANCERS): data/dictionary_long_names.json | $(OUTPUT_FOLDER)
	cp -f $< $@

INFO_DATASETS = $(patsubst %,%/info_datasets.csv, $(INPUTS))

$(STATS_COHORTS): $(INFO_DATASETS) $(PIPELINE)/cohorts.tsv data/dictionary_long_names.json | $(OUTPUT_FOLDER)
	python cohorts.py $@ $(PIPELINE)/cohorts.tsv data/dictionary_long_names.json $(INPUTS)



_GENES = $(BIOMART) $(CGC)
$(GENES)::  $(_GENES) | $(OUTPUT_FOLDER) $(OUTPUT_FOLDER_DRIVERS)
	python genes.py $@ $(BIOMART) $(CGC) \
		$(OUTPUT_FOLDER_DRIVERS)


$(DRIVER_GENES): $(PIPELINE)/drivers.tsv $(VEP) | $(OUTPUT_FOLDER)
	python drivers.py $@ $(PIPELINE)/drivers.tsv $(VEP)


$(MUTATIONS): $(PIPELINE)/mutations.tsv
	cat $< | gzip > $@


clean:
	@echo "Removing data"
	#rm -rf $(COUNT_SAMPLES) $(COUNT_VARIANTS)
	rm $(STATS_COHORTS) $(DRIVER_GENES) $(GENES) $(MUTATIONS)
	rm -rf $(OUTPUT_FOLDER_DRIVERS)
