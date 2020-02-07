SHELL := /bin/bash

GENOME = 38
ENSEMBL ?= 92
# TODO upgrade CADD version
CADD ?= 1.4
DATASETS?=../datasets/hg${GENOME}_ensembl${ENSEMBL}_$(shell date +%Y%m%d)
CORES ?= 1

# TODO
# export COSMIC_KEY=$(echo "email@example.com:mycosmicpassword" | base64)


# Configure BGData
BGDATA_LOCAL=$(DATASETS)/bgdata
BGDATA_OFFLINE="FALSE"
export BGDATA_LOCAL
export BGDATA_OFFLINE

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

# BGData package are managed by bgdata itself
.PHONY: bgdata
bgdata: | $(DATASETS)
	@echo bgdata datasets
	bgdata get datasets/genomereference/hg${GENOME}
	bgdata get intogen/expression/tcga_pancanatlas

ALL_TARGETS+=bgdata

####### SHARED #######
SHARED_FOLDER ?= $(DATASETS)/shared
$(SHARED_FOLDER):
	mkdir $@

TRANSCRIPTS_SQL_QUERY="SELECT g.stable_id, t.stable_id, x.display_label FROM gene g JOIN transcript t ON (g.canonical_transcript_id = t.transcript_id) JOIN xref x ON (g.display_xref_id = x.xref_id AND g.biotype='protein_coding') LEFT JOIN external_db ed USING (external_db_id) WHERE ed.db_name = 'HGNC';"

TRANSCRIPTS = $(SHARED_FOLDER)/ensembl_canonical_transcripts.tsv
$(TRANSCRIPTS): | $(SHARED_FOLDER)
	@echo Building ensembl canonical transcripts
	mysql -u anonymous -h ensembldb.ensembl.org --column-names=FALSE \
		-e ${TRANSCRIPTS_SQL_QUERY} ${ENSEMBL_DATABASE} > $@

ALL_TARGETS+=$(TRANSCRIPTS)

SHARED_BIOMART_QUERY=`cat shared/biomartQuery.txt`
SHARED_BIOMART_QUERY_ENCODED = $(shell python -c "from urllib.parse import quote_plus; query ='''${SHARED_BIOMART_QUERY}'''; print(quote_plus(query.replace('\n', '')))")
CDS_BIOMART = $(SHARED_FOLDER)/cds_biomart.tsv
$(CDS_BIOMART): $(TRANSCRIPTS) shared/biomartQuery.txt | $(SHARED_FOLDER)
	@echo Downloading biomart
	curl -s "${BIOMART_URL}?query=${SHARED_BIOMART_QUERY_ENCODED}" |\
		grep -f <(cut -f2 $(TRANSCRIPTS)) |\
		awk -F'\t' '($$5!=""){print($$0)}' > $@

ALL_TARGETS+=$(CDS_BIOMART)

REGIONS_CDS = $(SHARED_FOLDER)/cds.regions.gz
$(REGIONS_CDS): $(CDS_BIOMART) | $(SHARED_FOLDER)
	@echo Building CDS annotations
	echo -e "CHROMOSOME\tSTART\tEND\tSTRAND\tELEMENT\tSEGMENT\tSYMBOL" | \
		gzip > $@
	cat $(CDS_BIOMART) | \
		awk -F'\t' '($$5!=""){gsub("-1", "-", $$10); gsub("1", "+", $$10); print($$4"\t"$$5"\t"$$6"\t"$$10"\t"$$1"\t"$$1"\t"$$2)}' | \
		gzip >> $@

ALL_TARGETS+=$(REGIONS_CDS)

REGIONS_WG = $(SHARED_FOLDER)/wg.regions.gz
$(REGIONS_WG): $(SHARED_REGIONS_WG_SCRIPT) shared/create_wg_regions.py | $(SHARED_FOLDER)
	@echo Building whole-genome regions
	python shared/create_wg_regions.py hg${GENOME} 3 | gzip > $@

ALL_TARGETS+=$(REGIONS_WG)

COUNT_CDS = $(SHARED_FOLDER)/cds.counts.gz
$(COUNT_CDS): $(REGIONS_CDS) | $(SHARED_FOLDER)
	@echo Computing CDS signature
	bgsignature count -r $(REGIONS_CDS) -s 3 -g hg${GENOME} --cores ${CORES} --collapse --exclude-N -o $@

ALL_TARGETS+=$(COUNT_CDS)

COUNT_WG = $(SHARED_FOLDER)/wg.counts.gz
$(COUNT_WG): $(REGIONS_WG) | $(SHARED_FOLDER)
	@echo Computing whole-genome signature
	bgsignature count -r $(REGIONS_WG) -s 3 -g hg${GENOME} --cores ${CORES} --collapse --exclude-N -o $@

ALL_TARGETS+=$(COUNT_WG)

SOMATIC_PON_URL="https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd/download?path=%2FHMFTools-Resources%2FSage&files=SOMATIC_PON.vcf.gz"
SOMATIC_PON = $(SHARED_FOLDER)/somatic_pon_count_filtered.tsv.gz
$(SOMATIC_PON): shared/somatic_pon_counts.py | $(SHARED_FOLDER)
	@echo Getting somatic panel of normal counts
	python shared/somatic_pon_counts.py -u ${SOMATIC_PON_URL} -o $@

ALL_TARGETS+=$(SOMATIC_PON)

CONSEQUENCE_CDS = $(SHARED_FOLDER)/consequences.pickle.gz
# TODO check for other formats?
TRIPLETS_CDS = $(SHARED_FOLDER)/triplets.json.gz
# TODO this file exists in previous versions as triplets.pickle.gz
# TODO is it used?
# TODO avoid the use of bgvep
$(CONSEQUENCE_CDS) $(TRIPLETS_CDS) &: $(REGIONS_CDS) shared/count.py | $(SHARED_FOLDER)
	@echo Computing CDS triplets and consequences
	python shared/count.py -r $(REGIONS_CDS) -t ${TRIPLETS_CDS} -c $(CONSEQUENCE_CDS) \
		-v ${ENSEMBL} -g hg${GENOME}

ALL_TARGETS+=$(CONSEQUENCE_CDS)

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
GENOME_SIGNATURE = mutrate/signatures.cosmic.genome.tsv
MUTRATE_GENOME_SIGNATURE = $(MUTRATE_FOLDER)/signatures.cosmic.genome.tsv
$(MUTRATE_GENOME_SIGNATURE): $(GENOME_SIGNATURE) | $(MUTRATE_FOLDER)
	cp -f $< $@

ALL_TARGETS+=$(MUTRATE_GENOME_SIGNATURE)

EXOME_SIGNATURE = mutrate/signatures.cosmic.exome.tsv
MUTRATE_EXOME_SIGNATURE = $(MUTRATE_FOLDER)/signatures.cosmic.exome.tsv
$(MUTRATE_EXOME_SIGNATURE): $(EXOME_SIGNATURE) $(MUTRATE_GENOME_SIGNATURE) $(COUNT_CDS) $(COUNT_WG) mutrate/cosmic2exome.py | $(MUTRATE_FOLDER)
	@echo Building mutrate exome signature
	python mutrate/cosmic2exome.py $(EXOME_SIGNATURE) $(MUTRATE_GENOME_SIGNATURE) $(COUNT_CDS) $(COUNT_WG) $@

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


# TODO include FML configuration file


####### CBASE #######
CBASE_FOLDER = $(DATASETS)/cbase
$(CBASE_FOLDER): | $(DATASETS)
	mkdir $@

.PHONY: cbase

cbase: | $(CBASE_FOLDER)
	@echo Building CBaSE datasets
	wget -c http://genetics.bwh.harvard.edu/cbase/CBaSE_v1.1.zip -O $(CBASE_FOLDER)/cbase.zip
	unzip -d $(CBASE_FOLDER)/ $(CBASE_FOLDER)/cbase.zip
	mv $(CBASE_FOLDER)/CBaSE_v1.1/Auxiliary/*.gz $(CBASE_FOLDER)/
	mv $(CBASE_FOLDER)/CBaSE_v1.1/Auxiliary/*.txt $(CBASE_FOLDER)/
	rm -r $(CBASE_FOLDER)/CBaSE_v1.1
	rm -r $(CBASE_FOLDER)/__MACOSX/


ALL_TARGETS+=cbase

####### MUTPANNING #######
MUTPANNING_FOLDER = $(DATASETS)/mutpanning
$(MUTPANNING_FOLDER): | $(DATASETS)
	mkdir $@

.PHONY: mutpanning

mutpanning: | $(MUTPANNING_FOLDER)
	@echo Building MutPanning datasets
	wget https://datasets.genepattern.org/data/module_support_files/MutPanning/Hg19.zip -O $(MUTPANNING_FOLDER)/mutpanning.zip
	unzip -d $(MUTPANNING_FOLDER) $(MUTPANNING_FOLDER)/mutpanning.zip
	rm -r $(MUTPANNING_FOLDER)/__MACOSX/

ALL_TARGETS+=mutpanning

####### VEP #######
VEP_FOLDER = $(DATASETS)/vep
$(VEP_FOLDER): | $(DATASETS)
	mkdir $@

# TODO what is vep_install??


####### PTMS #######
# TODO not required?


####### BOOSTDM #######
# TODO not required?

####### COMBINATION #######
COMBINATION_FOLDER = $(DATASETS)/combination
$(COMBINATION_FOLDER): | $(DATASETS)
	mkdir $@

OLFACTORY_RECEPTORS = $(COMBINATION_FOLDER)/olfactory_receptors.tsv
$(OLFACTORY_RECEPTORS): | $(COMBINATION_FOLDER)
	wget https://genome.weizmann.ac.il/horde/download/genes.csv \
		-O $(OLFACTORY_RECEPTORS)

ALL_TARGETS += $(OLFACTORY_RECEPTORS)

COMBINATION_NEGATIVE_SET = $(COMBINATION_FOLDER)/negative_gene_set.tsv
COMBINATION_NON_EXPRESSED = $(COMBINATION_FOLDER)/non_expressed_genes_tcga.tsv
$(COMBINATION_NEGATIVE_SET) $(COMBINATION_NON_EXPRESSED) &: $(OLFACTORY_RECEPTORS) combination/create_negative_set.py | $(COMBINATION_FOLDER)
	@echo Building negative set
#	python combination/create_negative_set.py \
#		--olfactory_receptors $(OLFACTORY_RECEPTORS) \
#		--output_total $(COMBINATION_NEGATIVE_SET) \
#		--output_non_expressed $(COMBINATION_NON_EXPRESSED)

ALL_TARGETS += $(COMBINATION_NEGATIVE_SET) $(COMBINATION_NON_EXPRESSED)

GCC_FOLDER = $(COMBINATION_FOLDER)/cgc
$(GCC_FOLDER): | $(COMBINATION_FOLDER)
	mkdir $@

CGC = $(GCC_FOLDER)/cancer_gene_census.csv
$(CGC): combination/download_cgc.py | $(GCC_FOLDER)
	@echo Download CGC
	python combination/download_cgc.py --download $(GCC_FOLDER)

# TODO why 2 mappings
CGC_MAP=$(GCC_FOLDER)/mapping_cgc_ttypes.json
$(CGC_MAP): combination/mapping_cgc_ttypes.json
		cp -f $< $@

CGC_MAP_INTOGEN=$(GCC_FOLDER)/cgc/mapping_cgc_ttypes_intogen.json
$(CGC_MAP_INTOGEN): combination/mapping_cgc_ttypes_intogen.json
		cp -f $< $@

CGC_PARSED = $(GCC_FOLDER)/cancer_gene_census_parsed.csv
$(CGC_PARSED): $(CGC) $(CGC_MAP) $(CGC_MAP_INTOGEN) combination/parse_cgc.py | $(GCC_FOLDER)
	@echo Parsing CGC dataframe
#	python combination/parse_cgc.py \
#		--path_cgc_original $(CGC) \
#		--dict_mapping_cgc $(CGC_MAP) \
#		--dict_mapping_cgc_intogen $(CGC_MAP_INTOGEN) \
#		--path_output $(GCC_FOLDER) \

ALL_TARGETS += $(CGC) $(CGC_MAP) $(CGC_MAP_INTOGEN) $(CGC_PARSED)


#########################

all: $(ALL_TARGETS) | $(DATASETS)

clean:
	rm -rf $(DATASETS)
