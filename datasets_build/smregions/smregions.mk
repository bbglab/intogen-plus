
SRC_DATASETS_SMREGIONS = ${DATASETS_SOURCE_FOLDER}/smregions

DATASETS_SMREGIONS = $(DATASETS)/smregions
$(DATASETS_SMREGIONS): | $(DATASETS)
	mkdir $@


# Biomart Query
BIOMART_PFAM_QUERY=`cat ${SRC_DATASETS_SMREGIONS}/biomartQuery.txt`
BIOMART_PFAM_QUERY_ENCODED = $(shell python -c "from urllib.parse import quote_plus; query ='''${BIOMART_PFAM_QUERY}'''; print(quote_plus(query.replace('\n', '')))")
BIOMART_PFAM = $(DATASETS_SMREGIONS)/pfam_biomart.tsv.gz
$(BIOMART_PFAM): $$(TRANSCRIPTS) | $(DATASETS_SMREGIONS)
	@echo Downloading biomart
	curl -s "${BIOMART_URL}?query=${BIOMART_PFAM_QUERY_ENCODED}" |\
		grep -f <(cut -f2 $(TRANSCRIPTS)) |\
		awk -F'\t' '($$5!=""){print($$0)}' \
		| gzip > $@


REGIONS_PFAM = $(DATASETS_SMREGIONS)/regions_pfam.tsv
$(REGIONS_PFAM): $(BIOMART_PFAM) ${SRC_DATASETS_SMREGIONS}/panno.sh $$(CONTAINER_TRANSVAR) $$(TRANSCRIPTS) $$(DATASETS_TRANSVAR_FILES) | $(FOLDER)
	@echo Building CDS annotations
	echo -e "CHROMOSOME\tSTART\tEND\tSTRAND\tELEMENT_ID\tSEGMENT\tSYMBOL" \
		> $@
	zcat $(BIOMART_PFAM) | \
		awk '{system("${SRC_DATASETS_SMREGIONS}/panno.sh "$$2" "$$5" "$$3" "$$4" $(CONTAINER_TRANSVAR) $(DATASETS_TRANSVAR) $(TRANSCRIPTS)")}' \
		| grep -v "^\s" >> $@

TARGETS_DATASETS += $(BIOMART_PFAM) $(REGIONS_PFAM)