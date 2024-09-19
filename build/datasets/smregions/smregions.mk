
smregions_data_srcdir = ${src_datasets}/smregions

smregions_dir = $(INTOGEN_DATASETS)/smregions
$(smregions_dir): | $(INTOGEN_DATASETS)
	mkdir $@


# Biomart Query
biomart_pfram_query_file = ${smregions_data_srcdir}/biomartQuery.txt
biomart_pfram_query = `cat ${biomart_pfram_query_file}`
biomart_pfram_query_encoded = $(shell python -c "from urllib.parse import quote_plus; query ='''${biomart_pfram_query}'''; print(quote_plus(query.replace('\n', '')))")
BIOMART_PFAM = $(smregions_dir)/pfam_biomart.tsv.gz

$(BIOMART_PFAM): $(biomart_pfram_query_file) $$(ENSEMBL) | $(smregions_dir)
	@echo Downloading biomart
	@echo ${biomart_pfram_query}
	curl -L -s "${biomart_url}?query=${biomart_pfram_query_encoded}" |\
		tail -n +2 |\
		awk -F'\t' '($$4!=""){print($$0)}' \
		| gzip > $@


REGIONS_PFAM = $(smregions_dir)/regions_pfam.tsv
# TODO all transvar is requried for only this step. Can it be simplified?
$(REGIONS_PFAM): ${smregions_data_srcdir}/pfam.sh ${smregions_data_srcdir}/panno.sh $(BIOMART_PFAM) $$(TRANSVAR_CONTAINER) $$(TRANSVAR_FILES) $$(BIOMART_CDS) | $(smregions_dir)
	@echo Building CDS annotations
	$< $@ $(BIOMART_PFAM) ${smregions_data_srcdir}/panno.sh \
		$(TRANSVAR_CONTAINER) $(transvar_dir) $(BIOMART_CDS) ${cores}


DATASETS += $(BIOMART_PFAM) $(REGIONS_PFAM)