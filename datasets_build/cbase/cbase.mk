
DATASETS_CBASE = $(DATASETS)/cbase
CHECKPOINT_CBASE = ${DATASETS_CBASE}/.checkpoint
$(CHECKPOINT_CBASE):
	@echo Building CBaSE datasets
	mkdir -p ${DATASETS_CBASE}
	wget -c http://genetics.bwh.harvard.edu/cbase/CBaSE_v1.1.zip \
		-O ${DATASETS_CBASE}/cbase.zip
	unzip -d ${DATASETS_CBASE} ${DATASETS_CBASE}/cbase.zip
	mv ${DATASETS_CBASE}/CBaSE_v1.1/Auxiliary/*.gz ${DATASETS_CBASE}/
	mv ${DATASETS_CBASE}/CBaSE_v1.1/Auxiliary/*.txt ${DATASETS_CBASE}/
	rm -r ${DATASETS_CBASE}/__MACOSX/
	rm ${DATASETS_CBASE}/cbase.zip
	rm -r ${DATASETS_CBASE}/CBaSE_v1.1
	touch $@

TARGETS_DATASETS += $(CHECKPOINT_CBASE)