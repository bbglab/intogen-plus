
others_data_srcdir = ${src_datasets}/others

others_dir = $(INTOGEN_DATASETS)/others
$(others_dir): | $(INTOGEN_DATASETS)
	mkdir $@


somatic_pon_url = "https://nextcloud.hartwigmedicalfoundation.nl/s/LTiKTd8XxBqwaiC/download?path=%2FHMFTools-Resources%2FSage&files=SOMATIC_PON.vcf.gz"
SOMATIC_PON = $(others_dir)/somatic_pon_count_filtered.tsv.gz
$(SOMATIC_PON): ${others_data_srcdir}/somatic_pon_counts.py | $(others_dir)
	@echo Getting somatic panel of normal counts
	python $< -i ${somatic_pon_url} -o $@


OLFACTORY_RECEPTORS = $(others_dir)/olfactory_receptors.tsv

$(OLFACTORY_RECEPTORS): | $(others_dir)
	wget https://genome.weizmann.ac.il/horde/download/genes.csv \
		-O $(OLFACTORY_RECEPTORS)


NEGATIVE_GENE_SET = $(others_dir)/negative_gene_set.tsv
NON_EXPRESSED_GENES = $(others_dir)/non_expressed_genes_tcga.tsv

$(NEGATIVE_GENE_SET): ${others_data_srcdir}/create_negative_set.py $(OLFACTORY_RECEPTORS) | $(others_dir)
	@echo Building negative set
	python $< \
		--olfactory_receptors $(OLFACTORY_RECEPTORS) \
		--output_total $(NEGATIVE_GENE_SET) \
		--output_non_expressed $(NON_EXPRESSED_GENES)
	touch $(NEGATIVE_GENE_SET)
	touch $(NON_EXPRESSED_GENES)

$(NON_EXPRESSED_GENES): $(NEGATIVE_GENE_SET)
	# computed above
	$(NOOP)


DATASETS += $(SOMATIC_PON) $(OLFACTORY_RECEPTORS) \
	$(NEGATIVE_GENE_SET) $(NON_EXPRESSED_GENES)