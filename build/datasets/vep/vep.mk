
vep_data_srcdir = ${src_datasets}/vep

vep_dir = $(INTOGEN_DATASETS)/vep
$(vep_dir): | $(INTOGEN_DATASETS)
	mkdir $@


VEP_CACHE = $(vep_dir)/.vep${ensembl}_cache

$(VEP_CACHE): $$(VEP_CONTAINER) $$(GENOME) | $(vep_dir)
	@echo Building VEP datasets
	singularity exec -B $(vep_dir):/opt/vep/.vep $(VEP_CONTAINER) \
		perl /opt/vep/src/ensembl-vep/INSTALL.pl -a cf -s homo_sapiens -y ${grch} --c /opt/vep/.vep --CONVERT
	touch $@


VEP_MUTATIONS = $(vep_dir)/vep.tsv.gz
$(VEP_MUTATIONS): ${vep_data_srcdir}/run.sh $(vep_data_srcdir)/mutations.py $$(BOOST_ALL_REGION) $$(VEP_CONTAINER) $(VEP_CACHE) | $(vep_dir)
	$< $(BOOST_ALL_REGION) $(VEP_CONTAINER) $(vep_dir) $@ ${cores}


VEP_MUTATIONS_INDEX = $(vep_dir)/vep.tsv.gz.tbi
$(VEP_MUTATIONS_INDEX): $(VEP_MUTATIONS)
	tabix -f -s 1 -b 2 -e 2 $(VEP_MUTATIONS)

DATASETS += $(VEP_CACHE) $(VEP_MUTATIONS) $(VEP_MUTATIONS_INDEX)
