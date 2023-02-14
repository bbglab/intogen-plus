
others_data_srcdir = ${src_datasets}/others

others_dir = $(INTOGEN_DATASETS)/others
$(others_dir): | $(INTOGEN_DATASETS)
	mkdir $@


somatic_pon_url = "https://storage.googleapis.com/hmf-public/HMFtools-Resources/dna_pipeline/v5_31/38/variants/SageGermlinePon.98x.38.tsv.gz"
SOMATIC_PON = $(others_dir)/somatic_pon_count_filtered.tsv.gz
$(SOMATIC_PON): ${others_data_srcdir}/somatic_pon_counts.py | $(others_dir)
	@echo Getting somatic panel of normal counts
	python $< -i ${somatic_pon_url} -o $@

# # This is a temporal hack because the original file has been deleted # Update 20/01/21 msanchezg: Link has been restablished. The PON file has been updated. 
# SOMATIC_PON = $(others_dir)/somatic_pon_count_filtered.tsv.gz
# $(SOMATIC_PON): ${others_data_srcdir}/bgdata_copy.sh | $(others_dir)
	# $< $(others_dir)


OLFACTORY_RECEPTORS = $(others_dir)/olfactory_receptors.tsv

$(OLFACTORY_RECEPTORS): | $(others_dir)
	wget https://genome.weizmann.ac.il/horde/download/genes.csv \
		-O $(OLFACTORY_RECEPTORS)


NEGATIVE_GENE_SET = $(others_dir)/negative_gene_set.tsv
NON_EXPRESSED_GENES = $(others_dir)/non_expressed_genes_tcga.tsv

#Add mapping files
TTYPES_MAP = $(others_data_srcdir)/mapping_oncotree_ttypes.json
SYMBOLS_MAP = $(others_data_srcdir)/mapping_new_hugo_symbols.json


$(NEGATIVE_GENE_SET): ${others_data_srcdir}/create_negative_set.py $(OLFACTORY_RECEPTORS) | $(others_dir)
	@echo Building negative set
	python $< \
		--olfactory_receptors $(OLFACTORY_RECEPTORS) \
		--output_total $(NEGATIVE_GENE_SET) \
		--output_non_expressed $(NON_EXPRESSED_GENES) \
		--dict_mapping_ttypes $(TTYPES_MAP) \
		--dict_mapping_symbols $(SYMBOLS_MAP)
	touch $(NEGATIVE_GENE_SET)
	touch $(NON_EXPRESSED_GENES)

$(NON_EXPRESSED_GENES): $(NEGATIVE_GENE_SET)
	# computed above
	$(NOOP)


DATASETS += $(SOMATIC_PON) $(OLFACTORY_RECEPTORS) \
	$(NEGATIVE_GENE_SET) $(NON_EXPRESSED_GENES)
