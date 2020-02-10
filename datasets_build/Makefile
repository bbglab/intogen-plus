SHELL := /bin/bash

GENOME = 38
ENSEMBL ?= 92
# TODO upgrade CADD version
CADD ?= 1.4
DATASETS?=../datasets/hg${GENOME}_ensembl${ENSEMBL}_$(shell date +%Y%m%d)
CORES ?= 1

# TODO
# export COSMIC_KEY=$(echo "email@example.com:mycosmicpassword" | base64)


.PHONY: fall all clean

# accumulate all targets
ALL_TARGETS=

# useful variables
ENSEMBL_DATABASE="homo_sapiens_core_${ENSEMBL}_${GENOME}"
ENSEMBL_ARCHIVE = $(shell awk -v release=${ENSEMBL} '{if ($$1 == release) {print $$2}}' ensembl_archive.txt)
BIOMART_URL = http://${ENSEMBL_ARCHIVE}.archive.ensembl.org/biomart/martservice



# trick to make all the first target
# and set it to the end once all variables are defined
fall: all

# Ensure DATA_DIR exists before adding files
$(DATASETS):
	mkdir $@


####### bgdata #######
# bgdata is configured in this Makefile and not in its own,
# because this way the variables are exported to all datasets

# Configure BGData
BGDATA_LOCAL=$(DATASETS)/bgdata
BGDATA_OFFLINE="FALSE"
export BGDATA_LOCAL
export BGDATA_OFFLINE

# BGData is set as prerrequisite of the other SELF_CONTAINED targets
# to avoid issues when parallelizing the Makefile
.PHONY: bgdata
bgdata: | $(DATASETS)
	$(MAKE) -C $@ FOLDER=$(DATASETS)/$@ GENOME=${GENOME} ENSEMBL=${ENSEMBL}

ALL_TARGETS += bgdata

####### self-contained targets #######
SELF_CONTAINED = cbase mutpanning combination mutpanning shared
.PHONY: $(SELF_CONTAINED)

$(SELF_CONTAINED): bgdata | $(DATASETS)
	$(MAKE) -C $@ FOLDER=$(DATASETS)/$@ GENOME=${GENOME} ENSEMBL=${ENSEMBL}


ALL_TARGETS += $(SELF_CONTAINED)

# FIXME change the path to cgc to move it to the self-contained list
.PHONY: cgc
cgc: | $(DATASETS)
	$(MAKE) -C $@ FOLDER=$(DATASETS)/combination/$@

ALL_TARGETS += cgc



####### PREPROCESS #######
# TODO

####### TRANSVAR #######
TRANSVAR_FOLDER ?= $(DATASETS)/transvar
$(TRANSVAR_FOLDER):
	mkdir $@

# FIXME for other genomes than hg38, the GRCh version may differ
GENOME_FASTA = $(TRANSVAR_FOLDER)/Homo_sapiens.GRCh${GENOME}.fa
$(GENOME_FASTA): transvar/build_fasta.sh | $(TRANSVAR_FOLDER)
	@echo Build genome fasta file
	bash transvar/build_fasta.sh hg${GENOME} $@

ALL_TARGETS+=$(GENOME_FASTA)

ENSEMBL_GTF = $(TRANSVAR_FOLDER)/Homo_sapiens.GRCh${GENOME}.ENSEMBL.gtf.gz
$(ENSEMBL_GTF): | $(TRANSVAR_FOLDER)
	@echo Downloading ENSEMBL GTF
	wget "ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/homo_sapiens/${ENSEMBL_GTF}" \
		-O $@

ALL_TARGETS+=$(ENSEMBL_GTF)

# TODO index fasta and ensembl and configure ensemble. Requies transvar image


####### SMREGIONS #######
SMREGIONS_FOLDER ?= $(DATASETS)/smregions
$(SMREGIONS_FOLDER):
	mkdir $@

SMREGIONS_BIOMART_QUERY=`cat smregions/biomartQuery.txt`
SMREGIONS_BIOMART_QUERY_ENCODED = $(shell python -c "from urllib.parse import quote_plus; query ='''${SMREGIONS_BIOMART_QUERY}'''; print(quote_plus(query.replace('\n', '')))")
BIOMART_PFAM = $(SMREGIONS_FOLDER)/pfam_biomart.tsv
$(BIOMART_PFAM): $(TRANSCRIPTS) | $(SMREGIONS_FOLDER)
	@echo Downloading biomart
	curl -s "${BIOMART_URL}?query=${BIOMART_QUERY_ENCODED}" |\
		grep -f <(cut -f2 $(TRANSCRIPTS)) |\
		awk -F'\t' '($$5!=""){print($$0)}' > $@

ALL_TARGETS+=$(BIOMART_PFAM)

# TODO this guy needs transvar image
REGIONS_PFAM = $(SMREGIONS_FOLDER)/regions_pfam.tsv.gz

####### DNDSCV #######
DNDSCV_FOLDER ?= $(DATASETS)/dndscv
$(DNDSCV_FOLDER):
	mkdir $@

REF_RDA = $(DNDSCV_FOLDER)/RefCDS.rda
$(REF_RDA): $(GENOME_FASTA) | $(DNDSCV_FOLDER)
	@echo Building dNdSCV reference

# TODO this guy requires the dnds image


####### MUTRATE #######
MUTRATE_FOLDER ?= $(DATASETS)/mutrate
$(MUTRATE_FOLDER):
	mkdir $@

# TODO where are the signature files coming from
.PHONY: mutrate
mutrate:
	$(MAKE) -C mutrate FOLDER=$(DATASETS)/mutrate

ALL_TARGETS+=mutrate

MUTRATE_GENOME_SIGNATURE = $(MUTRATE_FOLDER)/signatures.cosmic.genome.tsv
MUTRATE_EXOME_SIGNATURE = $(MUTRATE_FOLDER)/signatures.cosmic.exome.tsv
$(MUTRATE_EXOME_SIGNATURE): mutrate $(MUTRATE_GENOME_SIGNATURE) $(COUNT_CDS) $(COUNT_WG) mutrate/cosmic2exome.py | $(MUTRATE_FOLDER)
	@echo Building mutrate exome signature
	python mutrate/cosmic2exome.py $(MUTRATE_GENOME_SIGNATURE) $(COUNT_CDS) $(COUNT_WG) $@

# TODO remove the original exome signature file as is not used

ALL_TARGETS+=$(MUTRATE_EXOME_SIGNATURE)


####### ONCODRIVEFML #######
ONCODRIVEFML_FOLDER ?= $(DATASETS)/oncodrivefml
$(ONCODRIVEFML_FOLDER):
	mkdir $@

# FIXME for other genomes than hg38, the GRCh version may differ
CADD_URL = http://krishna.gs.washington.edu/download/CADD/v${CADD}/GRCh${GENOME}/whole_genome_SNVs.tsv.gz
CADD_SCORES = $(ONCODRIVEFML_FOLDER)/cadd.tsv.gz
CADD_SCORES_INDEX = $(ONCODRIVEFML_FOLDER)/cadd.tsv.gz.tbi
$(CADD_SCORES) $(CADD_SCORES_INDEX) &: $(REGIONS_CDS) | $(ONCODRIVEFML_FOLDER)
	@echo Building OncodriveFML datasets
#	zcat $(REGIONS_CDS) | tail -n +2 |\
#		awk -v cadd="${CADD_URL}" '{system("tabix "cadd" "$$1":"$$2"-"$$3)}' |\
#			gzip > ${CADD_SCORES}.tmp
#	zcat ${CADD_SCORES}.tmp |\
#		sort --parallel=${CORES} -S 4G -k1,1 -k2,2n |\
#		uniq | bgzip > $(CADD_SCORES)
#	rm ${CADD_SCORES}.tmp
#	tabix -s 1 -b 2 -e 2 $(CADD_SCORES)


# TODO set a oneliner

ALL_TARGETS+=$(CADD_SCORES) $(CADD_SCORES_INDEX)



####### VEP #######
VEP_FOLDER = $(DATASETS)/vep
$(VEP_FOLDER): | $(DATASETS)
	mkdir $@

# TODO what is vep_install??


####### PTMS #######
# TODO not required?


####### BOOSTDM #######
# TODO not required?






#########################

all: $(SELF_CONTAINED) | $(DATASETS)

clean:
	rm -rf $(DATASETS)
