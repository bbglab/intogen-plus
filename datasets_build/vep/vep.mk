

VEP_DIR = $(DATASETS)/vep
VEP_CACHE = $(VEP_DIR)/.vep${ensmebl}_cache
$(VEP_CACHE): $$(VEP_CONTAINER) $$(GENOME)
	@echo Building VEP datasets
	mkdir -p $(VEP_DIR)
	singularity exec -B $(VEP_DIR):/opt/vep/.vep $(VEP_CONTAINER) \
		perl /opt/vep/src/ensembl-vep/INSTALL.pl -a cf -s homo_sapiens -y ${grch} --c /opt/vep/.vep --CONVERT
	touch $@

ALL_DATASETS += $(VEP_CACHE)