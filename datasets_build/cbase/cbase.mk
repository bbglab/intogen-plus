
FOLDER_CBASE = $(DATASETS)/cbase
CHECKPOINT_CBASE = ${FOLDER_CBASE}/.checkpoint
$(CHECKPOINT_CBASE):
	@echo Building CBaSE datasets
	mkdir -p ${FOLDER_CBASE}
	wget -c http://genetics.bwh.harvard.edu/cbase/CBaSE_v1.1.zip \
		-O ${FOLDER_CBASE}/cbase.zip
	unzip -d ${FOLDER_CBASE} ${FOLDER_CBASE}/cbase.zip
	mv ${FOLDER_CBASE}/CBaSE_v1.1/Auxiliary/*.gz ${FOLDER_CBASE}/
	mv ${FOLDER_CBASE}/CBaSE_v1.1/Auxiliary/*.txt ${FOLDER_CBASE}/
	rm -r ${FOLDER_CBASE}/__MACOSX/
	rm ${FOLDER_CBASE}/cbase.zip
	touch $@

ALL_TARGETS += $(CHECKPOINT_CBASE)