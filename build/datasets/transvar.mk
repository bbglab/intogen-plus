

transvar_dir ?= $(INTOGEN_DATASETS)/transvar
$(transvar_dir): | $(INTOGEN_DATASETS)
	mkdir $@


TRANSVAR_FASTA = $(transvar_dir)/Homo_sapiens.${grch}.fa

$(TRANSVAR_FASTA): $$(GENOME_FASTA) | $(transvar_dir)
	@echo Copy genome fasta file to transvar dir
	cp $< $@


TRANSVAR_FASTA_INDEX = $(TRANSVAR_FASTA).fai

$(TRANSVAR_FASTA_INDEX): $(TRANSVAR_FASTA) $$(TRANSVAR_CONTAINER) $$(GENOME) | $(transvar_dir)
	@echo Indexing genome fasta file
	singularity run -B $(transvar_dir):/data $(TRANSVAR_CONTAINER) \
		index --reference /data/$(notdir $<)
	singularity run -B $(transvar_dir):/data $(TRANSVAR_CONTAINER) \
		config -k reference -v /data/Homo_sapiens.${grch}.fa \
		--refversion ${grch}

ENSEMBL_GTF = $(transvar_dir)/Homo_sapiens.${grch}.${ensembl}.gtf.gz

$(ENSEMBL_GTF): $$(GENOME) $$(ENSEMBL) | $(transvar_dir)
	@echo Downloading ENSEMBL GTF
	wget "ftp://ftp.ensembl.org/pub/release-${ensembl}/gtf/homo_sapiens/Homo_sapiens.${grch}.${ensembl}.gtf.gz" \
		-O $@
	touch $@

ENSEMBL_INDEX = $(ENSEMBL_GTF).transvardb
$(ENSEMBL_INDEX): $(ENSEMBL_GTF) $(TRANSVAR_FASTA_INDEX) $$(TRANSVAR_CONTAINER) | $(transvar_dir)
	@echo Configure genome reference
	singularity run -B $(transvar_dir):/data $(TRANSVAR_CONTAINER) \
		index --ensembl /data/$(notdir $(ENSEMBL_GTF))
	singularity run -B $(transvar_dir):/data $(TRANSVAR_CONTAINER) \
		config -k ensembl -v /data/$(notdir $(ENSEMBL_INDEX)) \
		--refversion ${grch}


TRANSVAR_FILES = $(TRANSVAR_FASTA) $(TRANSVAR_FASTA_INDEX) $(ENSEMBL_GTF) $(ENSEMBL_INDEX)
DATASETS += $(TRANSVAR_FILES)