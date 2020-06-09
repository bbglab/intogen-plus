
dndscv_data_dir = $(INTOGEN_DATASETS)/dndscv

$(dndscv_data_dir): | $(INTOGEN_DATASETS)
	mkdir $@


REF_RDA = $(dndscv_data_dir)/RefCDS.rda

$(REF_RDA): $$(BIOMART_CDS) $$(GENOME_FASTA) $$(DNDSCV_CONTAINER) | $(dndscv_data_dir)
	@echo Building dNdSCV reference
	echo "library(dndscv); buildref(\"$(BIOMART_CDS)\", \"$(GENOME_FASTA)\", outfile = \"$@\")" | \
		singularity exec $(DNDSCV_CONTAINER) R --no-save


DATASETS += $(REF_RDA)