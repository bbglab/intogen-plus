
cgc_datasets_srcdir = ${src_datasets}/cgc

cgc_dir = $(INTOGEN_DATASETS)/cgc

$(cgc_dir): | $(INTOGEN_DATASETS)
	mkdir $@


CGC = $(cgc_dir)/cancer_gene_census.csv
$(CGC): ${cgc_datasets_srcdir}/download.py | $(cgc_dir)
	@echo Download CGC
ifeq ($(COSMIC_KEY), )
	$(error COSMIC_KEY not set)
endif
	python $< --download $(cgc_dir)


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



DATASETS += $(CGC) $(CGC_MAP) $(CGC_MAP_INTOGEN) $(CGC_PARSED)