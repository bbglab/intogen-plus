# TODO make this .mk files shell runable and set the DATASETS_SOURCE_FOLDER to be .. by default and some other trick
# like hotmaps-all or hotmaps-clean

tmpdir := $(shell mktemp -d)

SRC_DATASETS_HOTMAPS = ${DATASETS_SOURCE_FOLDER}/hotmaps

DATASETS_HOTMAPS = $(DATASETS)/hotmaps
$(DATASETS_HOTMAPS): | $(DATASETS)
	mkdir $@


HOTMAPS_PDB_INFO = $(DATASETS_HOTMAPS)/pdb_info.txt.gz
$(HOTMAPS_PDB_INFO): | $(DATASETS_HOTMAPS)
	wget http://karchinlab.org/data/HotMAPS/pdb_info.txt.gz -O $@

# FIXME this file is probably not needed and we only need sqlite DB
HOTMAPS_DB_DUMP = $(DATASETS_HOTMAPS)/mupit_modbase.sql
$(HOTMAPS_DB_DUMP): | $(DATASETS_HOTMAPS)
	wget http://karchinlab.org/data/HotMAPS/mupit_modbase.sql.gz -O $@.gz
	gunzip -c $@.gz > $@
	rm $@.gz

HOTMAPS_DB = $(DATASETS_HOTMAPS)/mupit_database.db
$(HOTMAPS_DB): $(HOTMAPS_DB_DUMP) ${SRC_DATASETS_HOTMAPS}/mysql2sqlite | $(DATASETS_HOTMAPS)
	${SRC_DATASETS_HOTMAPS}/mysql2sqlite $< | sqlite3 $@

## TODO maybe this one does not need a checkpoint because the rsync will do the work
#HOTMAPS_BIOUNIT = $(DATASETS_HOTMAPS)/.pdbbiounit.checkpoint
#$(HOTMAPS_BIOUNIT): $(HOTMAPS_PDB_INFO) | $(DATASETS_HOTMAPS)
#	mkdir -p $(DATASETS_HOTMAPS)/pdb/biounit/coordinates/all
#	mkdir -p ${tmpdir}/bioinfo
#	#wget -r ftp://ftp.wwpdb.org/pub/pdb/data/biounit/coordinates/all/ -P $(DATASETS_HOTMAPS)/pdb/biounit/coordinates/all
#	rsync -rptz --copy-links --port=33444 ftp.rcsb.org::ftp_data/biounit/coordinates/all/*.pdb1.gz ${tmpdir}/bioinfo
#	for it in `cut -f 1 $(HOTMAPS_PDB_INFO) | sort |uniq`; do
#		if [ -f $$it.pdb1.gz ]; then \
#			mv $$it.pdb1.gz $(DATASETS_HOTMAPS)/pdb/biounit/coordinates/all/
#		fi;
#	done
#	# rm -r ${tmpdir}/bioinfo
#	touch $@

## Extract only the used files instead of everything and then remove
#HOTMAPS_REFSEQ = $(DATASETS_HOTMAPS)/.pdbrefseq.checkpoint
#$(HOTMAPS_REFSEQ): | $(DATASETS_HOTMAPS)
#	mkdir -p $(DATASETS_HOTMAPS)/pdb/
#	wget ftp://salilab.org/databases/modbase/projects/genomes/H_sapiens/2013/ModBase_H_sapiens_2013_refseq.tar.xz -P $(DATASETS_HOTMAPS)/pdb/
#	tar -C $(DATASETS_HOTMAPS)/pdb/ -xf $(DATASETS_HOTMAPS)/pdb/ModBase_H_sapiens_2013_refseq.tar.xz
#	rm $(DATASETS_HOTMAPS)/pdb/ModBase_H_sapiens_2013_refseq.tar.xz
#	rm $(DATASETS_HOTMAPS)/pdb/ModBase_H_sapiens_2013_refseq/H_sapiens_2013_refseq.summary.txt
#	rm -r $(DATASETS_HOTMAPS)/pdb/ModBase_H_sapiens_2013_refseq/alignments
#	touch $@


#HOTMAPS_MODELS = $(DATASETS_HOTMAPS)/.pdbmodels.checkpoint
#$(HOTMAPS_MODELS): $(HOTMAPS_PDB_INFO) $(SRC_DATASETS_HOTMAPS)/models.sh | $(DATASETS_HOTMAPS)
#	mkdir -p $(DATASETS_HOTMAPS)/pdb/ModBase_H_sapiens_2013_GRCh37.70.pep.all/models/model
#	mkdir -p ${tmpdir}/pdbmodels
#	wget ftp://salilab.org/databases/modbase/projects/genomes/H_sapiens/2013/ModBase_H_sapiens_2013_GRCh37.70.pep.all.tar.xz -P ${tmpdir}/pdbmodels
#	tar -C ${tmpdir}/pdbmodels -xf ${tmpdir}/pdbmodels/ModBase_H_sapiens_2013_GRCh37.70.pep.all.tar.xz ModBase_H_sapiens_2013_GRCh37.70.pep.all/models/model
#	cd ${tmpdir}/pdbmodels/ModBase_H_sapiens_2013_GRCh37.70.pep.all/models/model
#	$(SRC_DATASETS_HOTMAPS)/models.sh $(HOTMAPS_PDB_INFO)
#	#rm -r ${tmpdir}/pdbmodels
#	touch $@


HOTMAPS_BIOUNIT = $(DATASETS_HOTMAPS)/.pdbbiounit.checkpoint
$(HOTMAPS_BIOUNIT): $(HOTMAPS_PDB_INFO) $(SRC_DATASETS_HOTMAPS)/biounit.sh | $(DATASETS_HOTMAPS)
	mkdir -p $(DATASETS_HOTMAPS)/pdb/biounit/coordinates/all
	$(SRC_DATASETS_HOTMAPS)/biounit.sh $(HOTMAPS_PDB_INFO) $(DATASETS_HOTMAPS)/pdb/biounit/coordinates/all
	touch $@

# Extract only the used files instead of everything and then remove
HOTMAPS_REFSEQ = $(DATASETS_HOTMAPS)/.pdbrefseq.checkpoint
$(HOTMAPS_REFSEQ): $(HOTMAPS_PDB_INFO) $(SRC_DATASETS_HOTMAPS)/ref.sh | $(DATASETS_HOTMAPS)
	mkdir -p $(DATASETS_HOTMAPS)/pdb/ModBase_H_sapiens_2013_refseq/models/model
	$(SRC_DATASETS_HOTMAPS)/ref.sh $(HOTMAPS_PDB_INFO) $(DATASETS_HOTMAPS)/pdb/ModBase_H_sapiens_2013_refseq/models/model
	touch $@


HOTMAPS_MODELS = $(DATASETS_HOTMAPS)/.pdbmodels.checkpoint
$(HOTMAPS_MODELS): $(HOTMAPS_PDB_INFO) $(SRC_DATASETS_HOTMAPS)/models.sh | $(DATASETS_HOTMAPS)
	mkdir -p $(DATASETS_HOTMAPS)/pdb/ModBase_H_sapiens_2013_GRCh37.70.pep.all/models/model
	$(SRC_DATASETS_HOTMAPS)/models.sh $(HOTMAPS_PDB_INFO) $(DATASETS_HOTMAPS)/pdb/ModBase_H_sapiens_2013_GRCh37.70.pep.all/models/model
	touch $@

TARGETS_DATASETS += $(HOTMAPS_DB_DUMP) $(HOTMAPS_DB) $(HOTMAPS_BIOUNIT) $(HOTMAPS_REFSEQ) $(HOTMAPS_MODELS)
