
MUTPANNING_DIR = $(DATASETS)/mutpanning
MUTPANNING_DATA = $(MUTPANNING_DIR)/.mutpanning
$(MUTPANNING_DATA):
	@echo Building MutPanning datasets
	mkdir -p $(MUTPANNING_DIR)
	wget https://datasets.genepattern.org/data/module_support_files/MutPanning/Hg19.zip \
		-O ${tmpdir}/mutpanning.zip
	unzip -d ${tmpdir} ${tmpdir}/mutpanning.zip
	mv ${tmpdir}/Hg19 ${MUTPANNING_DIR}/
	chmod -R g+r ${MUTPANNING_DIR}/Hg19
	touch $@

ALL_DATASETS += $(MUTPANNING_DATA)