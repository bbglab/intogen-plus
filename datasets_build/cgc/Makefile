
FOLDER ?= cgc

CGC = $(FOLDER)/cancer_gene_census.csv
# TODO why 2 mappings
CGC_MAP=$(FOLDER)/mapping_cgc_ttypes.json
CGC_MAP_INTOGEN=$(FOLDER)/cgc/mapping_cgc_ttypes_intogen.json
CGC_PARSED = $(FOLDER)/cancer_gene_census_parsed.csv

.PHONY: all clean

all: $(CGC) $(CGC_MAP) $(CGC_MAP_INTOGEN) $(CGC_PARSED)

$(FOLDER):
	mkdir $@

# TODO fix
$(CGC): download.py | $(FOLDER) check
	@echo Download CGC
	python download.py --download $(FOLDER)

$(FOLDER)/%.json: %.json | $(FOLDER)
	cp -f $< $@


$(CGC_PARSED): $(CGC) $(CGC_MAP) $(CGC_MAP_INTOGEN) parse.py | $(FOLDER)
	@echo Parsing CGC dataframe
#	python combination/parse.py \
#		--path_cgc_original $(CGC) \
#		--dict_mapping_cgc $(CGC_MAP) \
#		--dict_mapping_cgc_intogen $(CGC_MAP_INTOGEN) \
#		--path_output $(FOLDER) \

clean:
	rm -f $(CGC) $(CGC_MAP) $(CGC_MAP_INTOGEN) $(CGC_PARSED)


check:
ifeq (${COSMIC_KEY}, )
	$(error COSMIC_KEY not set)
endif