
cbase_dir = $(INTOGEN_DATASETS)/cbase
$(cbase_dir): | $(INTOGEN_DATASETS)
	mkdir $@

cbase_data_srcdir = ${src_datasets}/cbase

CBASE_DATA = ${cbase_dir}/.checkpoint

$(CBASE_DATA): ${cbase_data_srcdir}/download.sh | $(cbase_dir)
	@echo Building CBaSE datasets
	$< ${cbase_dir}
	touch $@

DATASETS += $(CBASE_DATA)