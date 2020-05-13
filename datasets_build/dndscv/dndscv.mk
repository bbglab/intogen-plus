
dndscv_datasets_srcdir = ${src_datasets}/dndscv

dndscv_dir = $(DATASETS)/dndscv

$(dndscv_dir): | $(DATASETS)
	mkdir $@


REF_RDA = $(dndscv_dir)/RefCDS.rda

$(REF_RDA): $$(BIOMART_CDS) $$(GENOME_FASTA) $$(DNDSCV_CONTAINER) | $(dndscv_dir)
	@echo Building dNdSCV reference
	echo "library(dndscv); buildref(\"$(BIOMART_CDS)\", \"$(GENOME_FASTA)\", outfile = \"$@\")" | \
		singularity exec $(DNDSCV_CONTAINER) R --no-save


DATASETS_TARGETS += $(REF_RDA)