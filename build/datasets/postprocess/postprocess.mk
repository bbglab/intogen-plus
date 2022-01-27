
postprocess_datasets_srcdir = ${src_datasets}/postprocess

postprocess_dir = $(INTOGEN_DATASETS)/postprocess

$(postprocess_dir): | $(INTOGEN_DATASETS)
	mkdir $@


# EXACT = $(postprocess_dir)/constraint.txt.gz
# $(EXACT): | $(postprocess_dir)
# 	wget -O $@ https://storage.googleapis.com/gnomad-public/release/2.1/ht/constraint/constraint.txt.bgz

ARTIFACTS = $(postprocess_dir)/artifacts.json
$(ARTIFACTS): ${postprocess_datasets_srcdir}/artifacts.json | $(postprocess_dir)
	cp $< $@

CANCERMINE = $(postprocess_dir)/cancermine_sentences.tsv
$(CANCERMINE): | $(postprocess_dir)
	wget -O $@ https://zenodo.org/record/5764207/files/cancermine_sentences.tsv #updated link (release 7/12/21)
	#https://zenodo.org/record/2662509/files/cancermine_sentences.tsv previous link


GENES_BLACKLIST = $(postprocess_dir)/black_listed.txt
GENES_WHITELIST = $(postprocess_dir)/white_listed.txt
$(postprocess_dir)/%_listed.txt: ${postprocess_datasets_srcdir}/%_listed.txt | $(postprocess_dir)
	cp -f $< $@


#Exclude MoA. This is automatically calculated from the cancer_gene_census_parsed.tsv file in the postprocessing step.
# GENES_MOA = $(postprocess_dir)/gene_MoA.tsv
# $(GENES_MOA): ${postprocess_datasets_srcdir}/gene_MoA.tsv | $(postprocess_dir)
	# cp -f $< $@


DATASETS += $(EXACT) $(ARTIFACTS) $(CANCERMINE) \
	$(GENES_BLACKLIST) $(GENES_WHITELIST) 
	#$(GENES_MOA)
	
