# TODO make this .mk files shell runable and set the DATASETS_SOURCE_FOLDER to be .. by default and some other trick
# like hotmaps-all or hotmaps-clean

tmpdir := $(shell mktemp -d)

hotpmaps_datasets_srcdir = ${src_datasets}/hotmaps

hotmaps_dir = $(INTOGEN_DATASETS)/hotmaps
$(hotmaps_dir): | $(INTOGEN_DATASETS)
	mkdir $@

# TODO use this file in the code if possible
HOTMAPS_PDB_INFO = $(hotmaps_dir)/pdb_info.txt.gz
$(HOTMAPS_PDB_INFO): | $(hotmaps_dir)
	wget http://karchinlab.org/data/HotMAPS/pdb_info.txt.gz -O $@

# FIXME this file is probably not needed and we only need sqlite DB
HOTMAPS_DB_DUMP = $(hotmaps_dir)/mupit_modbase.sql
$(HOTMAPS_DB_DUMP): | $(hotmaps_dir)
	wget http://karchinlab.org/data/HotMAPS/mupit_modbase.sql.gz -O $@.gz
	gunzip -c $@.gz > $@
	rm $@.gz

HOTMAPS_DB = $(hotmaps_dir)/mupit_database.db
$(HOTMAPS_DB): $(HOTMAPS_DB_DUMP) ${hotpmaps_datasets_srcdir}/mysql2sqlite | $(hotmaps_dir)
	${hotpmaps_datasets_srcdir}/mysql2sqlite $< | sqlite3 $@

## TODO maybe this one does not need a checkpoint because the rsync will do the work
#HOTMAPS_BIOUNIT = $(hotmaps_dir)/.pdbbiounit.checkpoint
#$(HOTMAPS_BIOUNIT): $(HOTMAPS_PDB_INFO) | $(hotmaps_dir)
#	mkdir -p $(hotmaps_dir)/pdb/biounit/coordinates/all
#	mkdir -p ${tmpdir}/bioinfo
#	#wget -r ftp://ftp.wwpdb.org/pub/pdb/data/biounit/coordinates/all/ -P $(hotmaps_dir)/pdb/biounit/coordinates/all
#	rsync -rptz --copy-links --port=33444 ftp.rcsb.org::ftp_data/biounit/coordinates/all/*.pdb1.gz ${tmpdir}/bioinfo
#	for it in `cut -f 1 $(HOTMAPS_PDB_INFO) | sort |uniq`; do
#		if [ -f $$it.pdb1.gz ]; then \
#			mv $$it.pdb1.gz $(hotmaps_dir)/pdb/biounit/coordinates/all/
#		fi;
#	done
#	# rm -r ${tmpdir}/bioinfo
#	touch $@

## Extract only the used files instead of everything and then remove
#HOTMAPS_REFSEQ = $(hotmaps_dir)/.pdbrefseq.checkpoint
#$(HOTMAPS_REFSEQ): | $(hotmaps_dir)
#	mkdir -p $(hotmaps_dir)/pdb/
#	wget ftp://salilab.org/databases/modbase/projects/genomes/H_sapiens/2013/ModBase_H_sapiens_2013_refseq.tar.xz -P $(hotmaps_dir)/pdb/
#	tar -C $(hotmaps_dir)/pdb/ -xf $(hotmaps_dir)/pdb/ModBase_H_sapiens_2013_refseq.tar.xz
#	rm $(hotmaps_dir)/pdb/ModBase_H_sapiens_2013_refseq.tar.xz
#	rm $(hotmaps_dir)/pdb/ModBase_H_sapiens_2013_refseq/H_sapiens_2013_refseq.summary.txt
#	rm -r $(hotmaps_dir)/pdb/ModBase_H_sapiens_2013_refseq/alignments
#	touch $@


#HOTMAPS_MODELS = $(hotmaps_dir)/.pdbmodels.checkpoint
#$(HOTMAPS_MODELS): $(HOTMAPS_PDB_INFO) $(hotpmaps_datasets_srcdir)/models.sh | $(hotmaps_dir)
#	mkdir -p $(hotmaps_dir)/pdb/ModBase_H_sapiens_2013_GRCh37.70.pep.all/models/model
#	mkdir -p ${tmpdir}/pdbmodels
#	wget ftp://salilab.org/databases/modbase/projects/genomes/H_sapiens/2013/ModBase_H_sapiens_2013_GRCh37.70.pep.all.tar.xz -P ${tmpdir}/pdbmodels
#	tar -C ${tmpdir}/pdbmodels -xf ${tmpdir}/pdbmodels/ModBase_H_sapiens_2013_GRCh37.70.pep.all.tar.xz ModBase_H_sapiens_2013_GRCh37.70.pep.all/models/model
#	cd ${tmpdir}/pdbmodels/ModBase_H_sapiens_2013_GRCh37.70.pep.all/models/model
#	$(hotpmaps_datasets_srcdir)/models.sh $(HOTMAPS_PDB_INFO)
#	#rm -r ${tmpdir}/pdbmodels
#	touch $@


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
	zcat $(HOTMAPS_PDB_INFO) | tail -n+2 | cut -f1 |sort |uniq > ${tmpdir}/list_pdbs.txt
	python $< $(HOTMAPS_DB) ${tmpdir}/list_pdbs.txt $@
	rm ${tmpdir}/list_pdbs.txt


# TODO used?
HOTMAPS_SIGNATURES = $(hotmaps_dir)/.signatures
$(HOTMAPS_SIGNATURES): $(hotpmaps_datasets_srcdir)/signatures.sh | $(hotmaps_dir)
	mkdir -p $(hotmaps_dir)/signatures
	$< $(hotmaps_dir)/signatures
	touch $@


# TODO get rid of this file
HOTMAPS_INFO_FULL = $(hotmaps_dir)/fully_described_pdb_info.txt
$(HOTMAPS_INFO_FULL): $(hotpmaps_datasets_srcdir)/info.sh | $(hotmaps_dir)
	$< $(hotmaps_dir)


# TODO get rid of this file
HOTMAPS_COORDINATES = $(hotmaps_dir)/coordinates.txt.gz
$(HOTMAPS_COORDINATES): $(hotpmaps_datasets_srcdir)/coordinates.sh | $(hotmaps_dir)
	$< $(hotmaps_dir)

DATASETS_TARGETS += $(HOTMAPS_DB_DUMP) $(HOTMAPS_DB) $(HOTMAPS_BIOUNIT) $(HOTMAPS_REFSEQ) $(HOTMAPS_MODELS) $(HOTMAPS_SIGNATURES) $(HOTMAPS_INFO_FULL) $(HOTMAPS_COORDINATES)
