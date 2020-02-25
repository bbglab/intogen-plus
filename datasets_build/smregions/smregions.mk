
SRC_DATASETS_SMREGIONS = ${DATASETS_SOURCE_FOLDER}/smregions

DATASETS_SMREGIONS = $(DATASETS)/smregions
$(DATASETS_SMREGIONS): | $(DATASETS)
	mkdir $@





# TODO requires transvar
REGIONS_PFAM = $(DATASETS_SMREGIONS)/regions_pfam.tsv.gz


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

#$(REGIONS_CDS): $(BIOMART) | $(FOLDER)
#	@echo Building CDS annotations
#	echo -e "CHROMOSOME\tSTART\tEND\tSTRAND\tELEMENT\tSEGMENT\tSYMBOL" | \
#		gzip > $@
#	cat $(BIOMART) | \
#		awk -F'\t' '($$5!=""){gsub("-1", "-", $$10); gsub("1", "+", $$10); print($$4"\t"$$5"\t"$$6"\t"$$10"\t"$$1"\t"$$1"\t"$$2)}' | \
#		gzip >> $@

TARGETS_DATASETS += $(BIOMART_PFAM)