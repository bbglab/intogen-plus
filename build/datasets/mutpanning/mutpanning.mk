
mutpanning_dir = $(INTOGEN_DATASETS)/mutpanning
$(mutpanning_dir): | $(INTOGEN_DATASETS)
	mkdir $@

mutpanning_data_srcdir = ${src_datasets}/mutpanning

MUTPANNING_DATA = ${mutpanning_dir}/.checkpoint

$(MUTPANNING_DATA): ${mutpanning_data_srcdir}/download.sh | $(mutpanning_dir)
	@echo Building MutPanning datasets
	$< ${mutpanning_dir}
	touch $@

DATASETS += $(MUTPANNING_DATA)
