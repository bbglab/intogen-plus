
FOLDER_CGC = $(DATASETS)/combination/cgc
$(FOLDER_CGC): | $(DATASETS)
	mkdir -p $@

CGC = $(FOLDER_CGC)/cancer_gene_census.csv
$(CGC): cgc/download.py | $(FOLDER_CGC) check-cosmic-key
	@echo Download CGC
	python cgc/download.py --download $(FOLDER_CGC)

# TODO why 2 mappings
CGC_MAP = $(FOLDER_CGC)/mapping_cgc_ttypes.json
CGC_MAP_INTOGEN = $(FOLDER_CGC)/mapping_cgc_ttypes_intogen.json
$(FOLDER_CGC)/%.json: cgc/%.json | $(FOLDER_CGC)
	cp -f $< $@

CGC_PARSED = $(FOLDER_CGC)/cancer_gene_census_parsed.tsv
$(CGC_PARSED): $(CGC) $(CGC_MAP) $(CGC_MAP_INTOGEN) cgc/parse.py | $(FOLDER_CGC)
	@echo Parsing CGC dataframe
	python cgc/parse.py \
		--path_cgc_original $(CGC) \
		--dict_mapping_cgc $(CGC_MAP) \
		--dict_mapping_cgc_intogen $(CGC_MAP_INTOGEN) \
		--path_output $(FOLDER_CGC)

ALL_TARGETS += $(CGC) $(CGC_MAP) $(CGC_MAP_INTOGEN) $(CGC_PARSED)

check-cosmic-key:
ifeq ($(COSMIC_KEY), )
	$(error COSMIC_KEY not set)
endif