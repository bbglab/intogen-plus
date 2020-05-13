
# Configure BGData
# Note that this configuration is applied to all datasets
BGDATA_LOCAL=$(DATASETS)/bgdata
BGDATA_OFFLINE="FALSE"
export BGDATA_LOCAL
export BGDATA_OFFLINE

.PHONY: bgdata
bgdata: $$(GENOME) | $(DATASETS)
	@echo Downloading bgdata datasets
	bgdata get datasets/genomereference/hg${genome}
	bgdata get intogen/expression/tcga_pancanatlas
	bgdata get intogen/coverage/hg${genome}

DATASETS_TARGETS += bgdata