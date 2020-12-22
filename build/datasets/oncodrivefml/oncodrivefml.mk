
fml_data_srcdir = ${src_datasets}/oncodrivefml


fml_dir = $(INTOGEN_DATASETS)/oncodrivefml
$(fml_dir): | $(INTOGEN_DATASETS)
	mkdir $@


CADD_URL = http://krishna.gs.washington.edu/download/CADD/v${CADD}/GRCh${GENOME}/whole_genome_SNVs.tsv.gz
#CADD_URL = /workspace/datasets/CADD/v${cadd}/hg${genome}/whole_genome_SNVs.tsv.gz
FML_SCORES = $(fml_dir)/scores.tsv.gz
$(FML_SCORES): ${fml_data_srcdir}/cadd.sh $$(REGIONS_CDS) $$(CADD) $$(GENOME) | $(fml_dir)
	@echo Building OncodriveFML datasets
	$< $(REGIONS_CDS) $(CADD_URL) ${cores} $@


FML_SCORES_INDEX = $(FML_SCORES).tbi
$(FML_SCORES_INDEX): $(FML_SCORES) | $(fml_dir)
	tabix -f -s 1 -b 2 -e 2 $(FML_SCORES)


DATASETS += $(FML_SCORES) $(FML_SCORES_INDEX)