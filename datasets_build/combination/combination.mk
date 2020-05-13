
combination_datasets_srcdir = ${src_datasets}/combination

combination_dir = $(DATASETS)/combination

$(combination_dir): | $(DATASETS)
	mkdir $@


OLFACTORY_RECEPTORS = $(combination_dir)/olfactory_receptors.tsv

$(OLFACTORY_RECEPTORS): | $(combination_dir)
	wget https://genome.weizmann.ac.il/horde/download/genes.csv \
		-O $(OLFACTORY_RECEPTORS)


NEGATIVE_GENE_SET = $(combination_dir)/negative_gene_set.tsv
NON_EXPRESSED_GENES = $(combination_dir)/non_expressed_genes_tcga.tsv

$(NEGATIVE_GENE_SET): ${combination_datasets_srcdir}/create_negative_set.py $(OLFACTORY_RECEPTORS) | $(combination_dir)
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


DATASETS_TARGETS += $(OLFACTORY_RECEPTORS) $(NEGATIVE_GENE_SET) $(NON_EXPRESSED_GENES)