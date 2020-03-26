
SRC_DATASETS_TRANSVAR = ${DATASETS_SOURCE_FOLDER}/transvar

DATASETS_TRANSVAR ?= $(DATASETS)/transvar
$(DATASETS_TRANSVAR): | $(DATASETS)
	mkdir $@

# FIXME for other genomes than hg38, the GRCh version may differ
GENOME_FASTA = $(DATASETS_TRANSVAR)/Homo_sapiens.GRCh${GENOME}.fa
$(GENOME_FASTA): ${SRC_DATASETS_TRANSVAR}/build_fasta.sh | $(DATASETS_TRANSVAR)
	@echo Build genome fasta file
	bash ${SRC_DATASETS_TRANSVAR}/build_fasta.sh hg${GENOME} $@


GENOME_FASTA_INDEX = $(GENOME_FASTA).fai
$(GENOME_FASTA_INDEX): $(GENOME_FASTA) $$(CONTAINER_TRANSVAR) | $(DATASETS_TRANSVAR)
	@echo Indexing genome fasta file
	singularity run -B $(DATASETS_TRANSVAR):/data $(CONTAINER_TRANSVAR) \
		index --reference /data/$$(basename $<)
	singularity run -B $(DATASETS_TRANSVAR):/data $(CONTAINER_TRANSVAR) \
		config -k reference -v /data/Homo_sapiens.GRCh${GENOME}.fa \
		--refversion GRCh${GENOME}

ENSEMBL_GTF = $(DATASETS_TRANSVAR)/Homo_sapiens.GRCh${GENOME}.${ENSEMBL}.gtf.gz
$(ENSEMBL_GTF): | $(DATASETS_TRANSVAR)
	@echo Downloading ENSEMBL GTF
	wget "ftp://ftp.ensembl.org/pub/release-${ENSEMBL}/gtf/homo_sapiens/Homo_sapiens.GRCh${GENOME}.${ENSEMBL}.gtf.gz" \
		-O $@

ENSEMBL_INDEX = $(ENSEMBL_GTF).transvardb
$(ENSEMBL_INDEX): $(ENSEMBL_GTF) $(GENOME_FASTA_INDEX) $$(CONTAINER_TRANSVAR) | $(DATASETS_TRANSVAR)
	@echo Configure genome reference
	singularity run -B $(DATASETS_TRANSVAR):/data $(CONTAINER_TRANSVAR) \
		index --ensembl /data/$$(basename $(ENSEMBL_GTF))
	singularity run -B $(DATASETS_TRANSVAR):/data $(CONTAINER_TRANSVAR) \
		config -k ensembl -v /data/$$(basename $(ENSEMBL_INDEX)) \
		--refversion GRCh${GENOME}

DATASETS_TRANSVAR_FILES = $(GENOME_FASTA) $(GENOME_FASTA_INDEX) $(ENSEMBL_GTF) $(ENSEMBL_INDEX)
TARGETS_DATASETS += $(GENOME_FASTA) $(GENOME_FASTA_INDEX) $(ENSEMBL_GTF) $(ENSEMBL_INDEX)