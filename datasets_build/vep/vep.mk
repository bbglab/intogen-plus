

DATASETS_VEP = $(DATASETS)/vep
CHECKPOINT_VEP = ${DATASETS_VEP}/.vep${ENSEMBL}.checkpoint
$(CHECKPOINT_VEP): $$(CONTAINER_VEP)
	@echo Building VEP datasets
	mkdir -p $(DATASETS_VEP)
	singularity exec -B $(DATASETS_VEP):/opt/vep/.vep $(CONTAINER_VEP) \
		perl /opt/vep/src/ensembl-vep/INSTALL.pl -a cf -s homo_sapiens -y GRCh${GENOME} --c /opt/vep/.vep --CONVERT
	touch $@

TARGETS_DATASETS += $(CHECKPOINT_VEP)