
FOLDER ?= mutrate

.PHONY: all clean

GENOME_SIGNATURE = $(FOLDER)/signatures.cosmic.genome.tsv

all: $(GENOME_SIGNATURE)

$(FOLDER):
	mkdir $@

$(GENOME_SIGNATURE): signatures.cosmic.genome.tsv | $(FOLDER)
	cp -f $< $@

clean:
	rm -f $(GENOME_SIGNATURE)
