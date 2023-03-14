
boostdm_dir = $(INTOGEN_DATASETS)/boostdm
$(boostdm_dir): | $$(INTOGEN_DATASETS)
	@echo Create directory
	mkdir -p $@

boostdm_data_src = $(src_datasets)/boostdm

BOOSTDM_PFAM = $(boostdm_dir)/pfam_biomart.tsv.gz

## STEP 2 - PFAM files
$(BOOSTDM_PFAM): $(boostdm_data_src)/symlink.sh $$(BIOMART_PFAM) $(boostdm_data_src)/download.sh | $(boostdm_dir) 
	@echo Creating biomart symlink
	$< $(boostdm_dir)/pfam_biomart.tsv.gz ${BIOMART_PFAM}
	$(boostdm_data_src)/download.sh $(boostdm_dir) 

## STEP 3 - ptms
BOOST_INFO_FUNCTIONAL = $(boostdm_dir)/ptms/info_functional_sites.json

$(BOOST_INFO_FUNCTIONAL): $(boostdm_data_src)/run_ptms.sh | $(boostdm_dir)
	@echo Creating phosphosite
	mkdir -p $(boostdm_dir)/ptms
	$< $(boostdm_dir)/ptms $(boostdm_data_src)

## STEP 4 - Symlinks
BOOST_SYMLINKS = $(boostdm_dir)/shared/.symlinks.checkpoint
SYMLINK_dir = $(boostdm_dir)/shared
$(BOOST_SYMLINKS): $(boostdm_data_src)/symlink.sh $$(REGIONS_CDS) $$(TRANSCRIPTS) $$(VEP_MUTATIONS) $$(VEP_MUTATIONS_INDEX) $$(oncotree_dir)| $(boostdm_dir)
	@echo Creating symlinks
	mkdir -p $(SYMLINK_dir)
	$^ $(SYMLINK_dir)
	touch $@
	
## STEP 4 - phylop 
BOOST_PHYLO = $(boostdm_dir)/hg38.phyloP100way.bw
$(BOOST_PHYLO): $(boostdm_data_src)/hg38.download.sh | $(boostdm_dir)
	@echo Creating phylop datasets
	$< $(boostdm_dir)

DATASETS += $(BOOSTDM_PFAM) $(BOOST_INFO_FUNCTIONAL) $(BOOST_PHYLO) $(BOOST_SYMLINKS)