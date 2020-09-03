# TODO where are the files coming from

# FIXME merge this with signature.mk into a single folder as this is used also by mutrate

deconstructsigs_data_srcdir = ${src_datasets}/deconstructsigs

deconstructsigs_dir = $(INTOGEN_DATASETS)/deconstructsigs

$(deconstructsigs_dir): | $(INTOGEN_DATASETS)
	mkdir $@


genome_signature = ${deconstructsigs_data_srcdir}/signatures.cosmic.genome.tsv


# TODO do not name it cosmic in the output
EXOME_SIGNATURE = $(deconstructsigs_dir)/signatures.cosmic.exome.tsv
$(EXOME_SIGNATURE): ${deconstructsigs_data_srcdir}/cosmic2exome.py $(MUTRATE_GENOME_SIGNATURE) $$(COUNT_CDS) $$(COUNT_WG) | $(deconstructsigs_dir)
	@echo Building mutrate exome signature
	python $< $(genome_signature) $(COUNT_CDS) $(COUNT_WG) $@

DATASETS += $(EXOME_SIGNATURE)