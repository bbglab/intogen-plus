
SRC_DATASETS_DNDSCV = ${DATASETS_SOURCE_FOLDER}/dndscv

DATASETS_DNDSCV = $(DATASETS)/dndscv
$(DATASETS_DNDSCV): | $(DATASETS)
	mkdir $@

REF_RDA = $(DATASETS_DNDSCV)/RefCDS.rda

$(REF_RDA): $$(BIOMART_CDS) $$(GENOME_FASTA) $$(CONTAINER_DNDSCV) | $(DATASETS_DNDSCV)
	@echo Building dNdSCV reference
	echo "library(dndscv); buildref(\"$(BIOMART_CDS)\", \"$(GENOME_FASTA)\", outfile = \"$@\")" | \
		singularity exec $(CONTAINER_DNDSCV) R --no-save


TARGETS_DATASETS += $(REF_RDA)