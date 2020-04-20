
smregions_data_srcdir = ${src_datasets}/smregions

SMREGIONS_DIR = $(DATASETS)/smregions
$(SMREGIONS_DIR): | $(DATASETS)
	mkdir $@


# Biomart Query
biomart_pfram_query_file = ${smregions_data_srcdir}/biomartQuery.txt
biomart_pfram_query = `cat ${biomart_pfram_query_file}`
biomart_pfram_query_encoded = $(shell python -c "from urllib.parse import quote_plus; query ='''${biomart_pfram_query}'''; print(quote_plus(query.replace('\n', '')))")
BIOMART_PFAM = $(SMREGIONS_DIR)/pfam_biomart.tsv.gz
$(BIOMART_PFAM): $$(TRANSCRIPTS) $(biomart_pfram_query_file) $$(ENSEMBL) | $(SMREGIONS_DIR)
	@echo Downloading biomart
	@echo ${biomart_pfram_query}
	curl -s "${biomart_url}?query=${biomart_pfram_query_encoded}" |\
		grep -f <(cut -f2 $(TRANSCRIPTS)) |\
		awk -F'\t' '($$5!=""){print($$0)}' \
		| gzip > $@


REGIONS_PFAM = $(SMREGIONS_DIR)/regions_pfam.tsv
$(REGIONS_PFAM): ${smregions_data_srcdir}/panno.sh $(BIOMART_PFAM) $$(CONTAINER_TRANSVAR) $$(TRANSCRIPTS) $$(TRANSVAR_FILES) | $(SMREGIONS_DIR)
	@echo Building CDS annotations
	echo -e "CHROMOSOME\tSTART\tEND\tSTRAND\tELEMENT_ID\tSEGMENT\tSYMBOL" \
		> $@
	zcat $(BIOMART_PFAM) | \
		awk '{system("$< "$$2" "$$5" "$$3" "$$4" $(CONTAINER_TRANSVAR) $(DATASETS_TRANSVAR) $(TRANSCRIPTS)")}' \
		| grep -v "^\s" >> $@

ALL_DATASETS += $(BIOMART_PFAM) $(REGIONS_PFAM)