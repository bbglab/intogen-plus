
FOLDER_ONCODRIVEFML = $(DATASETS)/oncodrivefml
$(FOLDER_ONCODRIVEFML):
	mkdir $@

# FIXME for other genomes than hg38, the GRCh version may differ
CADD_URL = http://krishna.gs.washington.edu/download/CADD/v${CADD}/GRCh${GENOME}/whole_genome_SNVs.tsv.gz
CADD_SCORES = ${FOLDER_ONCODRIVEFML}/cadd.tsv.gz
$(CADD_SCORES): $(REGIONS_CDS) | $(FOLDER_ONCODRIVEFML)
	@echo Building OncodriveFML datasets
#	mkdir -p ${ONCODRIVEFML}
#	zcat $(REGIONS_CDS) | tail -n +2 |\
#		awk -v cadd="${CADD_URL}" '{system("tabix "cadd" "$$1":"$$2"-"$$3)}' |\
#		gzip > ${CADD_SCORES}.tmp
#	zcat ${CADD_SCORES}.tmp |\
#		sort --parallel=${CORES} -S 4G -k1,1 -k2,2n |\
#		uniq | bgzip > $(CADD_SCORES)
#	rm ${CADD_SCORES}.tmp
# FIXME
	touch $@

# TODO set a oneliner

CADD_SCORES_INDEX = ${FOLDER_ONCODRIVEFML}/cadd.tsv.gz.tbi
$(CADD_SCORES_INDEX): | $(FOLDER_ONCODRIVEFML)
#	tabix -s 1 -b 2 -e 2 $(CADD_SCORES)
# FIXME
	touch $@

ALL_TARGETS += $(CADD_SCORES) $(CADD_SCORES_INDEX)