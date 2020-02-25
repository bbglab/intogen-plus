
SRC_DATASETS_COMBINATION = ${DATASETS_SOURCE_FOLDER}/combination

DATASETS_COMBINATION = $(DATASETS)/combination
$(DATASETS_COMBINATION): | $(DATASETS)
	mkdir $@

OLFACTORY_RECEPTORS = $(DATASETS_COMBINATION)/olfactory_receptors.tsv
$(OLFACTORY_RECEPTORS): | $(DATASETS_COMBINATION)
	wget https://genome.weizmann.ac.il/horde/download/genes.csv \
		-O $(OLFACTORY_RECEPTORS)

NEGATIVE_GENE_SET = $(DATASETS_COMBINATION)/negative_gene_set.tsv
NON_EXPRESSED_GENES = $(DATASETS_COMBINATION)/non_expressed_genes_tcga.tsv
$(NEGATIVE_GENE_SET): $(OLFACTORY_RECEPTORS) ${SRC_DATASETS_COMBINATION}/create_negative_set.py | $(FOLDER)
	@echo Building negative set
	python ${SRC_DATASETS_COMBINATION}/create_negative_set.py \
		--olfactory_receptors $(OLFACTORY_RECEPTORS) \
		--output_total $(NEGATIVE_GENE_SET) \
		--output_non_expressed $(NON_EXPRESSED_GENES)
	touch $(NEGATIVE_GENE_SET)
	touch $(NON_EXPRESSED_GENES)

$(NON_EXPRESSED_GENES): $(NEGATIVE_GENE_SET)
	# computed above
	$(NOOP)

TARGETS_DATASETS += $(OLFACTORY_RECEPTORS) $(NEGATIVE_GENE_SET) $(NON_EXPRESSED_GENES)