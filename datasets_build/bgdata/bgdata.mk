
# Configure BGData
BGDATA_LOCAL=$(DATASETS)/bgdata
BGDATA_OFFLINE="FALSE"
export BGDATA_LOCAL
export BGDATA_OFFLINE

.PHONY: bgdata
bgdata: | $(DATASETS)
	@echo Downloading bgdata datasets
	bgdata get datasets/genomereference/hg${GENOME}
	bgdata get intogen/expression/tcga_pancanatlas
	bgdata get intogen/coverage/hg${GENOME}

TARGETS_DATASETS += bgdata