.PHONY: all clean cohorts check-env

# The following environment variables must be defined
INTOGEN_DATASETS = ${INTOGEN_DATASETS_FOLDER}
DNDSCV_PAN_RUN = ${INTOGEN_DNDSCV_PAN_RUN}
INTOGEN_GENOME_VERSION = ${INTOGEN_GENOME}
INTOGEN_VEP_VERSION = ${INTOGEN_VEP}

CORES = 54

# FOLDERS with input datasets
INPUTS = ${INTOGEN_INPUTS}

# FOLDERS with the runs
RUNS = ${INTOGEN_RUNS}

EXPRESSION_TCGA = ${INTOGEN_DATASETS}/combination/non_expressed_genes_tcga.tsv
#ORFACTORY = ${INTOGEN_DATASETS}/combination/olfactory_receptors.tsv
# TODO fix. Should be in the one above
ORFACTORY = ${INTOGEN_DATASETS}/combination/olfactory_receptors.tsv
CGC = ${INTOGEN_DATASETS}/combination/cgc/cancer_gene_census_parsed.tsv
BIOMART = ${INTOGEN_DATASETS}/shared/cds_biomart.tsv
DNDSCV_PAN = ${DNDSCV_PAN_RUN}

# OUTPUTS
OUTPUT_FOLDER = output
CANCERS = $(OUTPUT_FOLDER)/ctypes.json
STATS_COHORTS = $(OUTPUT_FOLDER)/cohorts.tsv
GENES = $(OUTPUT_FOLDER)/genes.tsv
DRIVER_GENES = $(OUTPUT_FOLDER)/drivers.tsv
MUTATIONS = $(OUTPUT_FOLDER)/mutations.tsv.gz


# Other intermediate files
DRIVER_GENES_THRESHOLD = "05"
OUTPUT_FOLDER_DRIVERS = $(OUTPUT_FOLDER)/drivers
EXACT = $(OUTPUT_FOLDER_DRIVERS)/constraint.txt.gz
CANCERMINE = $(OUTPUT_FOLDER_DRIVERS)/cancermine_sentences.tsv

all: check-env cancers cohorts genes drivers mutations

cancers: $(CANCERS)
cohorts: check-env-datasets $(STATS_COHORTS)
genes: check-env-datasets $(GENES)
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

# TODO use count_samples.txt
COUNT_SAMPLES = $(patsubst %,%/count_samples.txt, $(INPUTS))
$(COUNT_SAMPLES): $(INFO_DATASETS)
	bgvariants groupby -q --cores $(CORES) -s "cut -f1 | uniq | sort -u | wc -l" -g DATASET $(dir $@) | awk '{print($$2" "$$1)}' | sort -n > $@

# TODO use count_variants.txt
COUNT_VARIANTS = $(patsubst %,%/count_variants.txt, $(INPUTS))
$(COUNT_VARIANTS): $(INFO_DATASETS)
	bgvariants count -q --cores $(CORES) -g DATASET $(dir $@) > $@


COHORT_FILES = $(COUNT_SAMPLES) $(COUNT_VARIANTS) $(INFO_DATASETS)
$(STATS_COHORTS): $(COHORT_FILES) | $(OUTPUT_FOLDER)
	python cohorts.py $@ $(INPUTS)



_GENES = $(BIOMART) $(CGC)
$(GENES)::  $(_GENES) | $(OUTPUT_FOLDER) $(OUTPUT_FOLDER_DRIVERS)
	python genes.py $@ $(BIOMART) $(CGC) \
		$(OUTPUT_FOLDER_DRIVERS)



$(EXACT): | $(OUTPUT_FOLDER_DRIVERS)
	wget -O $@ https://storage.googleapis.com/gnomad-public/release/2.1/ht/constraint/constraint.txt.bgz

$(CANCERMINE): | $(OUTPUT_FOLDER_DRIVERS)
	wget -O $@ https://zenodo.org/record/2662509/files/cancermine_sentences.tsv

DRIVER_GENES_FILES = $(wildcard $(patsubst %,%/deconstructsig/*.in.gz, $(RUNS)))
DRIVER_GENES_FILES += $(wildcard $(patsubst %,%/deconstructsig/*.signature_likelihood, $(RUNS)))
DRIVER_GENES_FILES += $(wildcard $(patsubst %,%/combination/*.$(DRIVER_GENES_THRESHOLD).out.gz, $(RUNS)))
DRIVER_GENES_FILES = $(wildcard $(patsubst %,%/smregions/*.out.gz, $(RUNS)))
DRIVER_GENES_FILES = $(wildcard $(patsubst %,%/oncodriveclustl/*.clusters.gz, $(RUNS)))
DRIVER_GENES_FILES = $(wildcard $(patsubst %,%/hotmaps/*.clusters.gz, $(RUNS)))
DRIVER_GENES_FILES = $(wildcard $(patsubst %,%/dndscv/*.out.gz, $(RUNS)))
DRIVER_GENES_FILES += $(EXPRESSION_TCGA) $(ORFACTORY) $(CANCERMINE) $(EXACT) $(CGC) $(ENSEMBL)
DRIVER_GENES_FILES += $(STATS_COHORTS) $(DNDSCV_PAN)
DRIVER_GENES_FILES += data/artifacts.json data/black_listed.txt data/white_listed.txt \
	data/gene_MoA.tsv
DRIVER_GENES_FILES += $(MUTATIONS)
$(DRIVER_GENES): $(DRIVER_GENES_FILES) | $(OUTPUT_FOLDER) $(OUTPUT_FOLDER_DRIVERS)
	python drivers.py $@ $(EXPRESSION_TCGA) $(ORFACTORY) $(CANCERMINE) \
		$(EXACT) $(CGC) $(BIOMART) $(STATS_COHORTS) $(DRIVER_GENES_THRESHOLD) \
		$(RUNS) $(OUTPUT_FOLDER_DRIVERS)

	python drivers_data.py $@ $(DNDSCV_PAN) $(INTOGEN_VEP_VERSION) $(INTOGEN_GENOME_VERSION) \
		$(MUTATIONS) $(RUNS) $(OUTPUT_FOLDER_DRIVERS)


MUTATION_FILES = $(wildcard $(patsubst %,%/vep/*.out.gz, $(RUNS)))
$(MUTATIONS): $(MUTATION_FILES) $(BIOMART)
	python mutations.py $@ $(BIOMART) $(RUNS)


clean:
	@echo "Removing data"
	#rm -rf $(COUNT_SAMPLES) $(COUNT_VARIANTS)
	rm $(STATS_COHORTS) $(DRIVER_GENES) $(GENES) $(MUTATIONS)
	rm -rf $(OUTPUT_FOLDER_DRIVERS)

check-env: check-env-datasets check-env-dndscv

check-env-datasets:
ifndef INTOGEN_DATASETS
	$(error INTOGEN_DATASETS is undefined)
endif

check-env-dndscv:
ifndef DNDSCV_PAN_RUN
	$(error DNDSCV_PAN_RUN is undefined)
endif

