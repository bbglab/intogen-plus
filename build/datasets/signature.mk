
signature_dir = $(INTOGEN_DATASETS)/signature
$(signature_dir): | $(INTOGEN_DATASETS)
	mkdir $@


COUNT_CDS = $(signature_dir)/cds.counts.gz

$(COUNT_CDS): $$(REGIONS_CDS) $$(GENOME) $$(SIGNATURE_CONTAINER) | $(signature_dir)
	@echo Computing CDS signature
	singularity exec $(SIGNATURE_CONTAINER) \
		bgsignature count -r $(REGIONS_CDS) -s 3 -g hg${genome} --cores ${cores} --collapse --exclude-N -o $@


COUNT_WG = $(signature_dir)/wg.counts.gz

$(COUNT_WG): $$(REGIONS_WG) $$(GENOME) $$(SIGNATURE_CONTAINER) | $(signature_dir)
	@echo Computing whole-genome signature
	singularity exec $(SIGNATURE_CONTAINER) \
		bgsignature count -r $(REGIONS_WG) -s 3 -g hg${genome} --cores ${cores} --collapse --exclude-N -o $@


DATASETS += $(COUNT_CDS) $(COUNT_WG)
