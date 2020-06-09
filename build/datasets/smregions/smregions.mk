
smregions_data_srcdir = ${src_datasets}/smregions

smregions_dir = $(INTOGEN_DATASETS)/smregions
$(smregions_dir): | $(INTOGEN_DATASETS)
	mkdir $@


# Biomart Query
biomart_pfram_query_file = ${smregions_data_srcdir}/biomartQuery.txt
biomart_pfram_query = `cat ${biomart_pfram_query_file}`
biomart_pfram_query_encoded = $(shell python -c "from urllib.parse import quote_plus; query ='''${biomart_pfram_query}'''; print(quote_plus(query.replace('\n', '')))")
BIOMART_PFAM = $(smregions_dir)/pfam_biomart.tsv.gz

$(BIOMART_PFAM): $$(TRANSCRIPTS) $(biomart_pfram_query_file) $$(ENSEMBL) | $(smregions_dir)
	@echo Downloading biomart
	@echo ${biomart_pfram_query}
	curl -s "${biomart_url}?query=${biomart_pfram_query_encoded}" |\
		grep -f <(cut -f2 $(TRANSCRIPTS)) |\
		awk -F'\t' '($$5!=""){print($$0)}' \
		| gzip > $@


REGIONS_PFAM = $(smregions_dir)/regions_pfam.tsv
# TODO all transvar is requried for only this step. Can it be simplified?
$(REGIONS_PFAM): ${smregions_data_srcdir}/panno.sh $(BIOMART_PFAM) $$(TRANSVAR_CONTAINER) $$(TRANSCRIPTS) $$(TRANSVAR_FILES) | $(smregions_dir)
	@echo Building CDS annotations
	echo -e "CHROMOSOME\tSTART\tEND\tSTRAND\tELEMENT_ID\tSEGMENT\tSYMBOL" \
		> $@
	zcat $(BIOMART_PFAM) | \
		awk '{system("$< "$$2" "$$5" "$$3" "$$4" $(TRANSVAR_CONTAINER) $(transvar_dir) $(TRANSCRIPTS)")}' \
		| grep -v "^\s" >> $@


DATASETS += $(BIOMART_PFAM) $(REGIONS_PFAM)