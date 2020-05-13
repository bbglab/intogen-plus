
fml_data_srcdir = ${src_datasets}/oncodrivefml


fml_dir = $(DATASETS)/oncodrivefml
$(fml_dir): | $(DATASETS)
	mkdir $@


# FIXME for other genomes than hg38, the GRCh version may differ
#CADD_URL = http://krishna.gs.washington.edu/download/CADD/v${CADD}/GRCh${GENOME}/whole_genome_SNVs.tsv.gz
# TODO fix this
CADD_URL = /workspace/datasets/CADD/v${cadd}/hg${genome}/whole_genome_SNVs.tsv.gz

CADD_SCORES = $(fml_dir)/cadd.tsv.gz
$(CADD_SCORES): ${fml_data_srcdir}/cadd.sh $$(REGIONS_CDS) $$(CADD)  $$(GENOME) | $(fml_dir)
	@echo Building OncodriveFML datasets
	$< $(REGIONS_CDS) $(CADD_URL) ${cores} $@


CADD_SCORES_INDEX = $(fml_dir)/cadd.tsv.gz.tbi

$(CADD_SCORES_INDEX): $(CADD_SCORES) | $(fml_dir)
	tabix -f -s 1 -b 2 -e 2 $(CADD_SCORES)


DATASETS_TARGETS += $(CADD_SCORES) $(CADD_SCORES_INDEX)