

FOLDER ?= bgdata

.PHONY: all

all: | check
	@echo Downloading bgdata datasets
	bgdata get datasets/genomereference/hg${GENOME}
	bgdata get intogen/expression/tcga_pancanatlas
	bgdata get intogen/coverage/hg${GENOME}


check:
ifeq ($(GENOME), )
	$(error GENOME not set)
endif