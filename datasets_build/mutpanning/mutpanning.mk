
mutpanning_dir = $(INTOGEN_DATASETS)/mutpanning

MUTPANNING_DATA = $(mutpanning_dir)/.mutpanning

$(MUTPANNING_DATA):
	@echo Building MutPanning datasets
	mkdir -p $(mutpanning_dir)
	wget https://datasets.genepattern.org/data/module_support_files/MutPanning/Hg19.zip \
		-O ${tmpdir}/mutpanning.zip
	unzip -d ${tmpdir} ${tmpdir}/mutpanning.zip
	mv ${tmpdir}/Hg19 ${mutpanning_dir}/
	chmod -R g+r ${mutpanning_dir}/Hg19
	touch $@


DATASETS_TARGETS += $(MUTPANNING_DATA)