
FOLDER_COMBINATION = $(DATASETS)/combination
$(FOLDER_COMBINATION): | $(DATASETS)
	mkdir $@

OLFACTORY_RECEPTORS = $(FOLDER_COMBINATION)/olfactory_receptors.tsv
$(OLFACTORY_RECEPTORS): | $(FOLDER_COMBINATION)
	wget https://genome.weizmann.ac.il/horde/download/genes.csv \
		-O $(OLFACTORY_RECEPTORS)

NEGATIVE_GENE_SET = $(FOLDER_COMBINATION)/negative_gene_set.tsv
NON_EXPRESSED_GENES = $(FOLDER_COMBINATION)/non_expressed_genes_tcga.tsv
$(NEGATIVE_GENE_SET) $(NON_EXPRESSED_GENES) &: $(OLFACTORY_RECEPTORS) combination/create_negative_set.py | $(FOLDER)
	@echo Building negative set
#	python combination/create_negative_set.py \
#		--olfactory_receptors $(OLFACTORY_RECEPTORS) \
#		--output_total $(NEGATIVE_GENE_SET) \
#		--output_non_expressed $(NON_EXPRESSED_GENES)

ALL_TARGETS += $(OLFACTORY_RECEPTORS) $(NEGATIVE_GENE_SET) $(NON_EXPRESSED_GENES)