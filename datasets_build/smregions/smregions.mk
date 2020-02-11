
FOLDER_SMREGIONS = $(DATASETS)/smregions
$(FOLDER_SMREGIONS): | $(DATASETS)
	mkdir $@





# TODO requires transvar
REGIONS_PFAM = $(FOLDER_SMREGIONS)/regions_pfam.tsv.gz


# Biomart Query
BIOMART_PFAM_QUERY=`cat smregions/biomartQuery.txt`
BIOMART_PFAM_QUERY_ENCODED = $(shell python -c "from urllib.parse import quote_plus; query ='''${BIOMART_PFAM_QUERY}'''; print(quote_plus(query.replace('\n', '')))")
BIOMART_PFAM = $(FOLDER_SMREGIONS)/pfam_biomart.tsv.gz
$(BIOMART_PFAM): | $(FOLDER_SMREGIONS)
	@echo Downloading biomart
	curl -s "${BIOMART_URL}?query=${BIOMART_PFAM_QUERY_ENCODED}" |\
		grep -f <(cut -f2 ${TRANSCRIPTS}) |\
		awk -F'\t' '($$5!=""){print($$0)}' \
		| gzip > $@

#$(REGIONS_CDS): $(BIOMART) | $(FOLDER)
#	@echo Building CDS annotations
#	echo -e "CHROMOSOME\tSTART\tEND\tSTRAND\tELEMENT\tSEGMENT\tSYMBOL" | \
#		gzip > $@
#	cat $(BIOMART) | \
#		awk -F'\t' '($$5!=""){gsub("-1", "-", $$10); gsub("1", "+", $$10); print($$4"\t"$$5"\t"$$6"\t"$$10"\t"$$1"\t"$$1"\t"$$2)}' | \
#		gzip >> $@

ALL_TARGETS += $(BIOMART_PFAM)