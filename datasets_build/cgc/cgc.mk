
cgc_datasets_srcdir = ${src_datasets}/cgc

cgc_dir = $(DATASETS)/combination/cgc
$(cgc_dir): | $(DATASETS)
	mkdir -p $@


CGC = $(cgc_dir)/cancer_gene_census.csv

$(CGC): ${cgc_datasets_srcdir}/download.py | $(cgc_dir) check-cosmic-key
	@echo Download CGC
	python $< --download $(cgc_dir)


# TODO why 2 mappings
CGC_MAP = $(cgc_dir)/mapping_cgc_ttypes.json
CGC_MAP_INTOGEN = $(cgc_dir)/mapping_cgc_ttypes_intogen.json

$(cgc_dir)/%.json: ${cgc_datasets_srcdir}/%.json | $(cgc_dir)
	cp -f $< $@


CGC_PARSED = $(cgc_dir)/cancer_gene_census_parsed.tsv

$(CGC_PARSED): ${cgc_datasets_srcdir}/parse.py $(CGC) $(CGC_MAP) $(CGC_MAP_INTOGEN) | $(cgc_dir)
	@echo Parsing CGC dataframe
	python $< \
		--path_cgc_original $(CGC) \
		--dict_mapping_cgc $(CGC_MAP) \
		--dict_mapping_cgc_intogen $(CGC_MAP_INTOGEN) \
		--path_output $(cgc_dir)


DATASETS_TARGETS += $(CGC) $(CGC_MAP) $(CGC_MAP_INTOGEN) $(CGC_PARSED)

check-cosmic-key:
	ifeq ($(COSMIC_KEY), )
		$(error COSMIC_KEY not set)
	endif