SHELL := /bin/bash

FOLDER ?= ${INTOGEN_DATASETS}/shared
ENSEMBL ?= ${INTOGEN_ENSEMBL}
GENOME ?= ${INTOGEN_GENOME}
CORES ?= 1

.PHONY: all clean

TRANSCRIPTS = $(FOLDER)/ensembl_canonical_transcripts.tsv
BIOMART = $(FOLDER)/cds_biomart.tsv
REGIONS_CDS = $(FOLDER)/cds.regions.gz
REGIONS_WG = $(FOLDER)/wg.regions.gz
COUNT_CDS = $(FOLDER)/cds.counts.gz
COUNT_WG = $(FOLDER)/wg.counts.gz
SOMATIC_PON = $(FOLDER)/somatic_pon_count_filtered.tsv.gz
CONSEQUENCE_CDS = $(FOLDER)/consequences.pickle.gz
# TODO check for other formats?
TRIPLETS_CDS = $(FOLDER)/triplets.json.gz
# TODO this file exists in previous versions as triplets.pickle.gz
# TODO is it used?


all: $(TRANSCRIPTS) $(BIOMART) $(REGIONS_CDS) $(REGIONS_WG) \
	$(COUNT_CDS) $(COUNT_WG) $(SOMATIC_PON) \
	$(CONSEQUENCE_CDS) $(TRIPLETS_CDS) | check


$(FOLDER):
	mkdir $@


# Ensembl transcripts
SQL_QUERY="SELECT g.stable_id, t.stable_id, x.display_label FROM gene g JOIN transcript t ON (g.canonical_transcript_id = t.transcript_id) JOIN xref x ON (g.display_xref_id = x.xref_id AND g.biotype='protein_coding') LEFT JOIN external_db ed USING (external_db_id) WHERE ed.db_name = 'HGNC';"
ENSEMBL_DATABASE="homo_sapiens_core_${ENSEMBL}_${GENOME}"


$(TRANSCRIPTS): | $(FOLDER)
	@echo Building ensembl canonical transcripts
	mysql -u anonymous -h ensembldb.ensembl.org --column-names=FALSE \
		-e ${SQL_QUERY} ${ENSEMBL_DATABASE} > $@


# Biomart file
ENSEMBL_ARCHIVE = $(shell awk -v release=${ENSEMBL} '{if ($$1 == release) {print $$2}}' ensembl_archive.txt)
BIOMART_URL = http://${ENSEMBL_ARCHIVE}.archive.ensembl.org/biomart/martservice
# Biomart Query
BIOMART_QUERY=`cat biomartQuery.txt`
BIOMART_QUERY_ENCODED = $(shell python -c "from urllib.parse import quote_plus; query ='''${BIOMART_QUERY}'''; print(quote_plus(query.replace('\n', '')))")

$(BIOMART): $(TRANSCRIPTS) | $(FOLDER)
	@echo Downloading biomart
	curl -s "${BIOMART_URL}?query=${BIOMART_QUERY_ENCODED}" |\
		grep -f <(cut -f2 $(TRANSCRIPTS)) |\
		awk -F'\t' '($$5!=""){print($$0)}' > $@

$(REGIONS_CDS): $(BIOMART) | $(FOLDER)
	@echo Building CDS annotations
	echo -e "CHROMOSOME\tSTART\tEND\tSTRAND\tELEMENT\tSEGMENT\tSYMBOL" | \
		gzip > $@
	cat $(BIOMART) | \
		awk -F'\t' '($$5!=""){gsub("-1", "-", $$10); gsub("1", "+", $$10); print($$4"\t"$$5"\t"$$6"\t"$$10"\t"$$1"\t"$$1"\t"$$2)}' | \
		gzip >> $@


$(REGIONS_WG): | $(FOLDER)
	@echo Building whole-genome regions
	python create_wg_regions.py hg${GENOME} 3 | gzip > $@

$(COUNT_CDS): $(REGIONS_CDS) | $(FOLDER)
	@echo Computing CDS signature
	bgsignature count -r $(REGIONS_CDS) -s 3 -g hg${INTOGEN_GENOME} --cores ${CORES} --collapse --exclude-N -o $@

$(COUNT_WG): $(REGIONS_WG) | $(FOLDER)
	@echo Computing whole-genome signature
	bgsignature count -r $(REGIONS_WG) -s 3 -g hg${INTOGEN_GENOME} --cores ${CORES} --collapse --exclude-N -o $@

SOMATIC_PON_URL="https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd/download?path=%2FHMFTools-Resources%2FSage&files=SOMATIC_PON.vcf.gz"

$(SOMATIC_PON): | $(FOLDER)
	@echo Getting somatic panel of normal counts
	python somatic_pon_counts.py -u ${SOMATIC_PON_URL} -o $@


# TODO avoid the use of bgvep
$(CONSEQUENCE_CDS) $(TRIPLETS_CDS) &: $(REGIONS_CDS) | $(FOLDER)
	@echo Computing CDS triplets and consequences
	python count.py -r $(REGIONS_CDS) -t ${TRIPLETS_CDS} -c $(CONSEQUENCE_CDS) \
		-v ${ENSEMBL} -g hg${GENOME}

clean:
	rm -rf $(FOLDER)

check:
ifeq ($(FOLDER), /shared)
	$(error FOLDER not set)
endif
ifeq (${ENSEMBL}, )
	$(error ENSEMBL not set)
endif
ifeq (${GENOME}, )
	$(error GENOME not set)
endif