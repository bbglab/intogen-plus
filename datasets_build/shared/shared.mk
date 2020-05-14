
shared_data_srcdir = ${src_datasets}/shared

shared_dir = $(INTOGEN_DATASETS)/shared
$(shared_dir): | $(INTOGEN_DATASETS)
	mkdir $@


# Ensembl transcripts
transcripts_sql_query = "SELECT g.stable_id, t.stable_id, x.display_label FROM gene g JOIN transcript t ON (g.canonical_transcript_id = t.transcript_id) JOIN xref x ON (g.display_xref_id = x.xref_id AND g.biotype='protein_coding') LEFT JOIN external_db ed USING (external_db_id) WHERE ed.db_name = 'HGNC';"
TRANSCRIPTS = $(shared_dir)/ensembl_canonical_transcripts.tsv

$(TRANSCRIPTS): $$(ENSEMBL) $$(GENOME) | $(shared_dir)
	@echo Building ensembl canonical transcripts
	mysql -u anonymous -h ensembldb.ensembl.org --column-names=FALSE \
		-e ${transcripts_sql_query} ${ensembl_db} > $@


# Biomart Query
biomart_cds_query_file = ${shared_data_srcdir}/biomartQuery.txt
biomart_cds_query = `cat ${biomart_cds_query_file}`
biomart_cds_query_encoded = $(shell python -c "from urllib.parse import quote_plus; query ='''${biomart_cds_query}'''; print(quote_plus(query.replace('\n', '')))")
BIOMART_CDS = $(shared_dir)/cds_biomart.tsv

$(BIOMART_CDS): $(TRANSCRIPTS) $(biomart_cds_query_file) $$(ENSEMBL) | $(shared_dir)
	@echo Downloading biomart
	curl -s "${biomart_url}?query=${biomart_cds_query_encoded}" |\
		grep -f <(cut -f2 $(TRANSCRIPTS)) |\
		awk -F'\t' '($$5!=""){print($$0)}' > $@


REGIONS_CDS = $(shared_dir)/cds.regions.gz

$(REGIONS_CDS): $(BIOMART_CDS) | $(shared_dir)
	@echo Building CDS annotations
	echo -e "CHROMOSOME\tSTART\tEND\tSTRAND\tELEMENT\tSEGMENT\tSYMBOL" | \
		gzip > $@
	cat $(BIOMART_CDS) | \
		awk -F'\t' '($$5!=""){gsub("-1", "-", $$10); gsub("1", "+", $$10); print($$4"\t"$$5"\t"$$6"\t"$$10"\t"$$1"\t"$$1"\t"$$2)}' | \
		gzip >> $@


REGIONS_WG = $(shared_dir)/wg.regions.gz

$(REGIONS_WG): ${shared_data_srcdir}/create_wg_regions.py $$(GENOME) | $(shared_dir)
	@echo Building whole-genome regions
	python $< hg${genome} 3 | gzip > $@


COUNT_CDS = $(shared_dir)/cds.counts.gz

$(COUNT_CDS): $(REGIONS_CDS) $$(GENOME) | $(shared_dir)
	@echo Computing CDS signature
	bgsignature count -r $(REGIONS_CDS) -s 3 -g hg${genome} --cores ${cores} --collapse --exclude-N -o $@


COUNT_WG = $(shared_dir)/wg.counts.gz

$(COUNT_WG): $(REGIONS_WG) $$(GENOME) | $(shared_dir)
	@echo Computing whole-genome signature
	bgsignature count -r $(REGIONS_WG) -s 3 -g hg${genome} --cores ${cores} --collapse --exclude-N -o $@


somatic_pon_url = "https://nextcloud.hartwigmedicalfoundation.nl/s/LTiKTd8XxBqwaiC/download?path=%2FHMFTools-Resources%2FSage&files=SOMATIC_PON.vcf.gz"
SOMATIC_PON = $(shared_dir)/somatic_pon_count_filtered.tsv.gz

$(SOMATIC_PON): ${shared_data_srcdir}/somatic_pon_counts.py | $(shared_dir)
	@echo Getting somatic panel of normal counts
	python $< -u ${somatic_pon_url} -o $@


CONSEQUENCE_CDS = $(shared_dir)/consequences.pickle.gz
# TODO check for other formats?
TRIPLETS_CDS = $(shared_dir)/triplets.json.gz
# TODO this file exists in previous versions as triplets.pickle.gz
# TODO is it used?
# TODO avoid the use of bgvep

$(CONSEQUENCE_CDS): ${shared_data_srcdir}/count.py $(REGIONS_CDS) $$(ENSEMBL) $$(GENOME) | $(shared_dir)
	@echo Computing CDS triplets and consequences
	python $< -r $(REGIONS_CDS) -t ${TRIPLETS_CDS} -c $(CONSEQUENCE_CDS) \
		-v ${ensembl} -g hg${genome}

# computed above
$(TRIPLETS_CDS): $(CONSEQUENCE_CDS)
	$(NOOP)

DATASETS_TARGETS += $(TRANSCRIPTS) $(REGIONS_CDS) $(REGIONS_WG) \
	$(COUNT_CDS) $(COUNT_WG) $(SOMATIC_PON) \
	$(CONSEQUENCE_CDS) $(TRIPLETS_CDS)
