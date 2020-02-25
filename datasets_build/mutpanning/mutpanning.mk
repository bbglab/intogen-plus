
DATASETS_MUTPANNING = $(DATASETS)/mutpanning
CHECKPOINT_MUTPANNING = ${DATASETS_MUTPANNING}/.checkpoint
$(CHECKPOINT_MUTPANNING):
	@echo Building MutPanning datasets
	mkdir -p ${DATASETS_MUTPANNING}
	wget https://datasets.genepattern.org/data/module_support_files/MutPanning/Hg19.zip \
		-O ${DATASETS_MUTPANNING}/mutpanning.zip
	unzip -d ${DATASETS_MUTPANNING} ${DATASETS_MUTPANNING}/mutpanning.zip
	rm -r ${DATASETS_MUTPANNING}/__MACOSX/
	rm ${DATASETS_MUTPANNING}/mutpanning.zip
	chmod -R g+r ${DATASETS_MUTPANNING}/Hg19
	touch $@

TARGETS_DATASETS += $(CHECKPOINT_MUTPANNING)