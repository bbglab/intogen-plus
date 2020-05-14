

vep_dir = $(INTOGEN_DATASETS)/vep
VEP_CACHE = $(vep_dir)/.vep${ensembl}_cache

$(VEP_CACHE): $$(VEP_CONTAINER) $$(GENOME)
	@echo Building VEP datasets
	mkdir -p $(vep_dir)
	singularity exec -B $(vep_dir):/opt/vep/.vep $(VEP_CONTAINER) \
		perl /opt/vep/src/ensembl-vep/INSTALL.pl -a cf -s homo_sapiens -y ${grch} --c /opt/vep/.vep --CONVERT
	touch $@


DATASETS_TARGETS += $(VEP_CACHE)