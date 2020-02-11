
FOLDER_MUTPANNING = $(DATASETS)/mutpanning
CHECKPOINT_MUTPANNING = ${FOLDER_MUTPANNING}/.checkpoint
$(CHECKPOINT_MUTPANNING):
	@echo Building MutPanning datasets
	mkdir -p ${FOLDER_MUTPANNING}
	wget https://datasets.genepattern.org/data/module_support_files/MutPanning/Hg19.zip \
		-O ${FOLDER_MUTPANNING}/mutpanning.zip
	unzip -d ${FOLDER_MUTPANNING} ${FOLDER_MUTPANNING}/mutpanning.zip
	rm -r ${FOLDER_MUTPANNING}/__MACOSX/
	rm ${FOLDER_MUTPANNING}/mutpanning.zip
	chmod -R g+r ${FOLDER_MUTPANNING}/Hg19
	touch $@

ALL_TARGETS += $(CHECKPOINT_MUTPANNING)