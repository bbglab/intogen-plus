
transvar_data_srcdir = ${src_datasets}/transvar

TRANSVAR_DIR ?= $(DATASETS)/transvar
$(TRANSVAR_DIR): | $(DATASETS)
	mkdir $@

# FIXME for other genomes than hg38, the GRCh version may differ
genome_fasta_script = ${transvar_data_srcdir}/build_fasta.sh
GENOME_FASTA = $(TRANSVAR_DIR)/Homo_sapiens.${grch}.fa
$(GENOME_FASTA): $(genome_fasta_script) $$(GENOME) | $(TRANSVAR_DIR)
	@echo Build genome fasta file
	bash $< hg${genome} $@


GENOME_FASTA_INDEX = $(GENOME_FASTA).fai
$(GENOME_FASTA_INDEX): $(GENOME_FASTA) $$(TRANSVAR_CONTAINER) | $(TRANSVAR_DIR)
	@echo Indexing genome fasta file
	singularity run -B $(TRANSVAR_DIR):/data $(TRANSVAR_CONTAINER) \
		index --reference /data/$(notdir $<)
	singularity run -B $(TRANSVAR_DIR):/data $(TRANSVAR_CONTAINER) \
		config -k reference -v /data/Homo_sapiens.${grch}.fa \
		--refversion ${grch}

ENSEMBL_GTF = $(TRANSVAR_DIR)/Homo_sapiens.${grch}.${ensmebl}.gtf.gz
$(ENSEMBL_GTF): $$(GENOME) $$(ENSEMBL) | $(TRANSVAR_DIR)
	@echo Downloading ENSEMBL GTF
	wget "ftp://ftp.ensembl.org/pub/release-${ensmebl}/gtf/homo_sapiens/Homo_sapiens.${grch}.${ensmebl}.gtf.gz" \
		-O $@
	touch $@

ENSEMBL_INDEX = $(ENSEMBL_GTF).transvardb
$(ENSEMBL_INDEX): $(ENSEMBL_GTF) $(GENOME_FASTA_INDEX) $$(TRANSVAR_CONTAINER) | $(TRANSVAR_DIR)
	@echo Configure genome reference
	singularity run -B $(TRANSVAR_DIR):/data $(TRANSVAR_CONTAINER) \
		index --ensembl /data/$(notdir $(ENSEMBL_GTF))
	singularity run -B $(TRANSVAR_DIR):/data $(TRANSVAR_CONTAINER) \
		config -k ensembl -v /data/$(notdir $(ENSEMBL_INDEX)) \
		--refversion ${grch}

TRANSVAR_FILES = $(GENOME_FASTA) $(GENOME_FASTA_INDEX) $(ENSEMBL_GTF) $(ENSEMBL_INDEX)
TARGETS_DATASETS += $(TRANSVAR_FILES)