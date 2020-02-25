
SRC_TRANSVAR = ${DATASETS_SOURCE_FOLDER}/transvar

TRANSVAR_FOLDER ?= $(DATASETS)/transvar
$(TRANSVAR_FOLDER): | $(DATASETS)
	mkdir $@

# FIXME for other genomes than hg38, the GRCh version may differ
GENOME_FASTA = $(TRANSVAR_FOLDER)/Homo_sapiens.GRCh${GENOME}.fa
$(GENOME_FASTA): ${SRC_TRANSVAR}/build_fasta.sh | $(TRANSVAR_FOLDER)
	@echo Build genome fasta file
	bash ${SRC_TRANSVAR}/build_fasta.sh hg${GENOME} $@



ENSEMBL_GTF = $(TRANSVAR_FOLDER)/Homo_sapiens.GRCh${GENOME}.ENSEMBL.gtf.gz
$(ENSEMBL_GTF): | $(TRANSVAR_FOLDER)
	@echo Downloading ENSEMBL GTF
	wget "ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/homo_sapiens/${ENSEMBL_GTF}" \
		-O $@

#ALL_TARGETS+=$(ENSEMBL_GTF)

# TODO index fasta and ensembl and configure ensemble. Requies transvar image

TARGETS_DATASETS += $(GENOME_FASTA)