
FOLDER ?= combination

.PHONY: all clean

OLFACTORY_RECEPTORS = $(FOLDER)/olfactory_receptors.tsv
NEGATIVE_SET = $(FOLDER)/negative_gene_set.tsv
NON_EXPRESSED = $(FOLDER)/non_expressed_genes_tcga.tsv

all: $(OLFACTORY_RECEPTORS) $(NEGATIVE_SET) $(NON_EXPRESSED)

$(FOLDER):
	mkdir $@

$(OLFACTORY_RECEPTORS): | $(FOLDER)
	wget https://genome.weizmann.ac.il/horde/download/genes.csv \
		-O $(OLFACTORY_RECEPTORS)


$(NEGATIVE_SET) $(NON_EXPRESSED) &: $(OLFACTORY_RECEPTORS) create_negative_set.py | $(FOLDER)
	@echo Building negative set
#	python create_negative_set.py \
#		--olfactory_receptors $(OLFACTORY_RECEPTORS) \
#		--output_total $(NEGATIVE_SET) \
#		--output_non_expressed $(NON_EXPRESSED)

clean:
	rm -f $(OLFACTORY_RECEPTORS) $(NEGATIVE_SET) $(NON_EXPRESSED)
