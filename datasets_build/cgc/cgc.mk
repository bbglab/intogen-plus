
SRC_CGC = ${DATASETS_SOURCE_FOLDER}/cgc
DATASETS_CGC = $(DATASETS)/combination/cgc
$(DATASETS_CGC): | $(DATASETS)
	mkdir -p $@

CGC = $(DATASETS_CGC)/cancer_gene_census.csv
$(CGC): ${SRC_CGC}/download.py | $(DATASETS_CGC) check-cosmic-key
	@echo Download CGC
	python ${SRC_CGC}/download.py --download $(DATASETS_CGC)

# TODO why 2 mappings
CGC_MAP = $(DATASETS_CGC)/mapping_cgc_ttypes.json
CGC_MAP_INTOGEN = $(DATASETS_CGC)/mapping_cgc_ttypes_intogen.json
$(DATASETS_CGC)/%.json: ${SRC_CGC}/%.json | $(DATASETS_CGC)
	cp -f $< $@

CGC_PARSED = $(DATASETS_CGC)/cancer_gene_census_parsed.tsv
$(CGC_PARSED): $(CGC) $(CGC_MAP) $(CGC_MAP_INTOGEN) ${SRC_CGC}/parse.py | $(DATASETS_CGC)
	@echo Parsing CGC dataframe
	python ${SRC_CGC}/parse.py \
		--path_cgc_original $(CGC) \
		--dict_mapping_cgc $(CGC_MAP) \
		--dict_mapping_cgc_intogen $(CGC_MAP_INTOGEN) \
		--path_output $(DATASETS_CGC)

TARGETS_DATASETS += $(CGC) $(CGC_MAP) $(CGC_MAP_INTOGEN) $(CGC_PARSED)

check-cosmic-key:
ifeq ($(COSMIC_KEY), )
	$(error COSMIC_KEY not set)
endif