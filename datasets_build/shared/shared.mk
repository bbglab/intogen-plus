
SRC_DATASETS_SHARED = ${DATASETS_SOURCE_FOLDER}/shared

DATASETS_SHARED = $(DATASETS)/shared
$(DATASETS_SHARED): | $(DATASETS)
	mkdir $@

# Ensembl transcripts
TRANSCRIPTS_SQL_QUERY="SELECT g.stable_id, t.stable_id, x.display_label FROM gene g JOIN transcript t ON (g.canonical_transcript_id = t.transcript_id) JOIN xref x ON (g.display_xref_id = x.xref_id AND g.biotype='protein_coding') LEFT JOIN external_db ed USING (external_db_id) WHERE ed.db_name = 'HGNC';"
TRANSCRIPTS = $(DATASETS_SHARED)/ensembl_canonical_transcripts.tsv
$(TRANSCRIPTS): | $(DATASETS_SHARED)
	@echo Building ensembl canonical transcripts
	mysql -u anonymous -h ensembldb.ensembl.org --column-names=FALSE \
		-e ${TRANSCRIPTS_SQL_QUERY} ${ENSEMBL_DATABASE} > $@


# Biomart Query
BIOMART_CDS_QUERY=`cat ${SRC_DATASETS_SHARED}/biomartQuery.txt`
BIOMART_CDS_QUERY_ENCODED = $(shell python -c "from urllib.parse import quote_plus; query ='''${BIOMART_CDS_QUERY}'''; print(quote_plus(query.replace('\n', '')))")
BIOMART_CDS = $(DATASETS_SHARED)/cds_biomart.tsv
$(BIOMART_CDS): $(TRANSCRIPTS) | $(DATASETS_SHARED)
	@echo Downloading biomart
	curl -s "${BIOMART_URL}?query=${BIOMART_CDS_QUERY_ENCODED}" |\
		grep -f <(cut -f2 $(TRANSCRIPTS)) |\
		awk -F'\t' '($$5!=""){print($$0)}' > $@

REGIONS_CDS = $(DATASETS_SHARED)/cds.regions.gz
$(REGIONS_CDS): $(BIOMART_CDS) | $(DATASETS_SHARED)
	@echo Building CDS annotations
	echo -e "CHROMOSOME\tSTART\tEND\tSTRAND\tELEMENT\tSEGMENT\tSYMBOL" | \
		gzip > $@
	cat $(BIOMART_CDS) | \
		awk -F'\t' '($$5!=""){gsub("-1", "-", $$10); gsub("1", "+", $$10); print($$4"\t"$$5"\t"$$6"\t"$$10"\t"$$1"\t"$$1"\t"$$2)}' | \
		gzip >> $@

REGIONS_WG = $(DATASETS_SHARED)/wg.regions.gz
$(REGIONS_WG): ${SRC_DATASETS_SHARED}/create_wg_regions.py | $(DATASETS_SHARED)
	@echo Building whole-genome regions
	python ${SRC_DATASETS_SHARED}/create_wg_regions.py hg${GENOME} 3 | gzip > $@

COUNT_CDS = $(DATASETS_SHARED)/cds.counts.gz
$(COUNT_CDS): $(REGIONS_CDS) | $(DATASETS_SHARED)
	@echo Computing CDS signature
	bgsignature count -r $(REGIONS_CDS) -s 3 -g hg${GENOME} --cores ${CORES} --collapse --exclude-N -o $@

COUNT_WG = $(DATASETS_SHARED)/wg.counts.gz
$(COUNT_WG): $(REGIONS_WG) | $(DATASETS_SHARED)
	@echo Computing whole-genome signature
	bgsignature count -r $(REGIONS_WG) -s 3 -g hg${GENOME} --cores ${CORES} --collapse --exclude-N -o $@

SOMATIC_PON_URL="https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd/download?path=%2FHMFTools-Resources%2FSage&files=SOMATIC_PON.vcf.gz"
SOMATIC_PON = $(DATASETS_SHARED)/somatic_pon_count_filtered.tsv.gz
$(SOMATIC_PON): ${SRC_DATASETS_SHARED}/somatic_pon_counts.py | $(DATASETS_SHARED)
	@echo Getting somatic panel of normal counts
	python ${SRC_DATASETS_SHARED}/somatic_pon_counts.py -u ${SOMATIC_PON_URL} -o $@

CONSEQUENCE_CDS = $(DATASETS_SHARED)/consequences.pickle.gz
# TODO check for other formats?
TRIPLETS_CDS = $(DATASETS_SHARED)/triplets.json.gz
# TODO this file exists in previous versions as triplets.pickle.gz
# TODO is it used?
# TODO avoid the use of bgvep
$(CONSEQUENCE_CDS): $(REGIONS_CDS) ${SRC_DATASETS_SHARED}/count.py | $(DATASETS_SHARED)
	@echo Computing CDS triplets and consequences
	python ${SRC_DATASETS_SHARED}/count.py -r $(REGIONS_CDS) -t ${TRIPLETS_CDS} -c $(CONSEQUENCE_CDS) \
		-v ${ENSEMBL} -g hg${GENOME}

# computed above
$(TRIPLETS_CDS): $(CONSEQUENCE_CDS)
	$(NOOP)

TARGETS_DATASETS += $(TRANSCRIPTS) $(BIOMART) $(REGIONS_CDS) $(REGIONS_WG) \
	$(COUNT_CDS) $(COUNT_WG) $(SOMATIC_PON) \
	$(CONSEQUENCE_CDS) $(TRIPLETS_CDS)
