
FOLDER ?= mutpanning

.PHONY: all clean

CHECKPOINT = ${FOLDER}/.checkpoint

all: $(CHECKPOINT)


$(CHECKPOINT):
	@echo Building MutPanning datasets
	mkdir -p ${FOLDER}
	wget https://datasets.genepattern.org/data/module_support_files/MutPanning/Hg19.zip \
		-O ${FOLDER}/mutpanning.zip
	unzip -d ${FOLDER} ${FOLDER}/mutpanning.zip
	rm -r ${FOLDER}/__MACOSX/
	rm ${FOLDER}/mutpanning.zip
	chmod -R g+r ${FOLDER}/Hg19
	touch $@


clean:
	rm -f $(CHECKPOINT)
