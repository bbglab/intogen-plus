
DATASETS_FML = $(DATASETS)/oncodrivefml
$(DATASETS_FML): | $(DATASETS)
	mkdir $@

# FIXME for other genomes than hg38, the GRCh version may differ
#CADD_URL = http://krishna.gs.washington.edu/download/CADD/v${CADD}/GRCh${GENOME}/whole_genome_SNVs.tsv.gz
# TODO fix this
CADD_URL = /workspace/datasets/CADD/v1.4/hg38/whole_genome_SNVs.tsv.gz
CADD_SCORES = ${DATASETS_FML}/cadd.tsv.gz
$(CADD_SCORES): $$(REGIONS_CDS) | $(DATASETS_FML)
	@echo Building OncodriveFML datasets
	zcat $(REGIONS_CDS) | tail -n +2 |\
		awk -v cadd="${CADD_URL}" '{system("tabix "cadd" "$$1":"$$2"-"$$3)}' |\
		awk 'BEGIN {FS="\t";OFS = FS};{ $$5=""; print }' |\
		gzip > $@.tmp
	zcat $@.tmp |\
		sort --parallel=${CORES} -S 4G -k1,1 -k2,2n |\
		uniq | bgzip > $@
	rm $@.tmp


# TODO set a oneliner

CADD_SCORES_INDEX = $(DATASETS_FML)/cadd.tsv.gz.tbi
$(CADD_SCORES_INDEX): $(CADD_SCORES) | $(DATASETS_FML)
	tabix -s 1 -b 2 -e 2 $(CADD_SCORES)

TARGETS_DATASETS += $(CADD_SCORES) $(CADD_SCORES_INDEX)