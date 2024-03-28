
regions_data_srcdir = ${src_datasets}/regions

regions_dir = $(INTOGEN_DATASETS)/regions
$(regions_dir): | $(INTOGEN_DATASETS)
	mkdir $@


# Biomart Query
biomart_cds_query_file = ${regions_data_srcdir}/biomartQuery.txt
biomart_cds_query = `cat ${biomart_cds_query_file}`
biomart_cds_query_encoded = $(shell python -c "from urllib.parse import quote_plus; query ='''${biomart_cds_query}'''; print(quote_plus(query.replace('\n', '')))")
BIOMART_CDS = $(regions_dir)/cds_biomart.tsv

$(BIOMART_CDS): $(biomart_cds_query_file) $$(ENSEMBL) | $(regions_dir)
	@echo Downloading biomart
	curl -L -s "${biomart_url}?query=${biomart_cds_query_encoded}" |\
		tail -n +2 |\
		awk -F'\t' '($$5!=""){print($$0)}' > $@


REGIONS_CDS = $(regions_dir)/cds.regions.gz

$(REGIONS_CDS): $(BIOMART_CDS) | $(regions_dir)
	@echo Building CDS annotations
	echo -e "CHROMOSOME\tSTART\tEND\tSTRAND\tELEMENT\tSEGMENT\tSYMBOL" | \
		gzip > $@
	cat $(BIOMART_CDS) | \
		awk -F'\t' '($$5!=""){gsub("-1", "-", $$10); gsub("1", "+", $$10); print($$4"\t"$$5"\t"$$6"\t"$$10"\t"$$1"\t"$$1"\t"$$2)}' | \
		gzip >> $@


REGIONS_WG = $(regions_dir)/wg.regions.gz

$(REGIONS_WG): ${regions_data_srcdir}/create_wg_regions.py $$(GENOME) | $(regions_dir)
	@echo Building whole-genome regions
	python $< hg${genome} 3 | gzip > $@


GENOME_FASTA = $(regions_dir)/Homo_sapiens.${grch}.fa
$(GENOME_FASTA): ${regions_data_srcdir}/build_fasta.sh $$(GENOME) | $(regions_dir)
	@echo Build genome fasta file
	bash $< hg${genome} $@

DATASETS += $(REGIONS_CDS) $(REGIONS_WG) \
	$(GENOME_FASTA)
