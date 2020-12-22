
# Configure BGData
# Note that this configuration is applied to all datasets
BGDATA_LOCAL=$(INTOGEN_DATASETS)/bgdata
BGDATA_OFFLINE="FALSE"
export BGDATA_LOCAL
export BGDATA_OFFLINE

.PHONY: bgdata
bgdata: | $(INTOGEN_DATASETS)
	@echo Downloading bgdata datasets
	bgdata get datasets/genomereference/hg19
	bgdata get datasets/genomereference/hg38
	bgdata get intogen/coverage/hg19
	bgdata get intogen/coverage/hg38
	bgdata get intogen/dndscv/pan
