
oncotree_datasets_srcdir = ${src_datasets}/oncotree

oncotree_dir = $(INTOGEN_DATASETS)/oncotree

$(oncotree_dir): | $(INTOGEN_DATASETS)
	mkdir $@


ONCOTREE = $(oncotree_dir)/tree.tsv
$(ONCOTREE): ${oncotree_datasets_srcdir}/tree.tsv | $(oncotree_dir)
	cp $< $@

ONCOTREE_DEF = $(oncotree_dir)/oncotree-definitions.json
$(ONCOTREE_DEF): $(oncotree_datasets_srcdir)/dictionary_long_names.json | $(oncotree_dir)
	cp $< $@

DATASETS += $(ONCOTREE) $(ONCOTREE_DEF)