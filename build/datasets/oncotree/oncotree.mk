
oncotree_datasets_srcdir = ${src_datasets}/oncotree

oncotree_dir = $(INTOGEN_DATASETS)/oncotree

$(oncotree_dir): | $(INTOGEN_DATASETS)
	mkdir $@


ONCOTREE = $(oncotree_dir)/tree.tsv
$(ONCOTREE): ${oncotree_datasets_srcdir}/tree.tsv | $(oncotree_dir)
	cp $< $@


DATASETS += $(ONCOTREE)