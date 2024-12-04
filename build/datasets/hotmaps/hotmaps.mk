
hotpmaps_datasets_srcdir = ${src_datasets}/hotmaps

hotmaps_dir = $(INTOGEN_DATASETS)/hotmaps
$(hotmaps_dir): | $(INTOGEN_DATASETS)
	mkdir $@

# TODO use this file in the code if possible
HOTMAPS_PDB_INFO = $(hotmaps_dir)/pdb_info.txt.gz
$(HOTMAPS_PDB_INFO): | $(hotmaps_dir)

	wget 'https://www.dropbox.com/scl/fi/0lklxk9h8fkzzwrz9ge1x/pdb_info.txt.gz?rlkey=uadoiylkv1pcuaed267q751eu&st=n1ds6yhq&dl=1' -O $@

HOTMAPS_DB = $(hotmaps_dir)/mupit_database.db
$(HOTMAPS_DB): $(HOTMAPS_DB_DUMP) ${hotpmaps_datasets_srcdir}/mysql2sqlite | $(hotmaps_dir)
	wget 'https://www.dropbox.com/scl/fi/2fvqta3mcpy0jjs8ta3fc/mupit_modbase.sql.gz?rlkey=z2ud0kaqyr1p23kryuxv9o168&st=21vhlb4f&dl=1' -O $@.mysql.gz
	gunzip -c $@.mysql.gz > $@.mysql
	${hotpmaps_datasets_srcdir}/mysql2sqlite $@.mysql | sqlite3 $@
	rm $@.mysql.gz
	rm $@.mysql


HOTMAPS_BIOUNIT = $(hotmaps_dir)/.pdbbiounit.checkpoint
$(HOTMAPS_BIOUNIT): $(HOTMAPS_PDB_INFO) $(hotpmaps_datasets_srcdir)/biounit.sh | $(hotmaps_dir)
	mkdir -p $(hotmaps_dir)/pdb/biounit/coordinates/all
	$(hotpmaps_datasets_srcdir)/biounit.sh $(HOTMAPS_PDB_INFO) $(hotmaps_dir)/pdb/biounit/coordinates/all
	touch $@

# Extract only the used files instead of everything and then remove
HOTMAPS_REFSEQ = $(hotmaps_dir)/.pdbrefseq.checkpoint
$(HOTMAPS_REFSEQ): $(HOTMAPS_PDB_INFO) $(hotpmaps_datasets_srcdir)/ref.sh | $(hotmaps_dir)
	mkdir -p $(hotmaps_dir)/pdb/ModBase_H_sapiens_2013_refseq/models/model
	$(hotpmaps_datasets_srcdir)/ref.sh $(HOTMAPS_PDB_INFO) $(hotmaps_dir)/pdb/ModBase_H_sapiens_2013_refseq/models/model
	touch $@


HOTMAPS_MODELS = $(hotmaps_dir)/.pdbmodels.checkpoint
$(HOTMAPS_MODELS): $(HOTMAPS_PDB_INFO) $(hotpmaps_datasets_srcdir)/models.sh | $(hotmaps_dir)
	mkdir -p $(hotmaps_dir)/pdb/ModBase_H_sapiens_2013_GRCh37.70.pep.all/models/model
	$(hotpmaps_datasets_srcdir)/models.sh $(HOTMAPS_PDB_INFO) $(hotmaps_dir)/pdb/ModBase_H_sapiens_2013_GRCh37.70.pep.all/models/model
	touch $@


HOTMAPS_COORDINATES = $(hotmaps_dir)/coordinates.txt
$(HOTMAPS_COORDINATES): $(hotpmaps_datasets_srcdir)/generate_coordinates.py $(HOTMAPS_PDB_INFO) $(HOTMAPS_DB) | $(hotmaps_dir)
	zcat $(HOTMAPS_PDB_INFO) | tail -n+2 | cut -f1 |sort |uniq > $(@D)/list_pdbs.txt
	python $< $(HOTMAPS_DB) $(@D)/list_pdbs.txt $@
	rm $(@D)/list_pdbs.txt


# TODO get rid of this file
HOTMAPS_INFO_FULL = $(hotmaps_dir)/fully_described_pdb_info.txt
$(HOTMAPS_INFO_FULL): $(hotpmaps_datasets_srcdir)/info.sh | $(hotmaps_dir)
	$< $(hotmaps_dir)


#HOTMAPS_COORDINATES = $(hotmaps_dir)/coordinates.txt.gz
#$(HOTMAPS_COORDINATES): $(hotpmaps_datasets_srcdir)/generate_coordinates.py $(HOTMAPS_PDB_INFO) $(HOTMAPS_DB) | $(hotmaps_dir)
#	zcat $(HOTMAPS_PDB_INFO) | tail -n+2 | cut -f1 |sort |uniq > ${tmpdir}/list_pdbs.txt
#	python $< $(HOTMAPS_DB) ${tmpdir}/list_pdbs.txt $@
#	rm ${tmpdir}/list_pdbs.txt

HOTMAPS_COORDINATES = $(hotmaps_dir)/coordinates.txt.gz
$(HOTMAPS_COORDINATES): $(hotpmaps_datasets_srcdir)/coordinates.sh | $(hotmaps_dir)
	$< $(hotmaps_dir)

DATASETS += $(HOTMAPS_DB) $(HOTMAPS_COORDINATES) \
	$(HOTMAPS_INFO_FULL) $(HOTMAPS_COORDINATES) \
	$(HOTMAPS_BIOUNIT) $(HOTMAPS_REFSEQ) $(HOTMAPS_MODELS) $(HOTMAPS_SIGNATURES)
