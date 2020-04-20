
DATASETS_CBASE = $(DATASETS)/cbase
CBASE_DATA = ${DATASETS_CBASE}/.checkpoint
$(CBASE_DATA):
	@echo Building CBaSE datasets
	mkdir -p ${DATASETS_CBASE}
	wget -c http://genetics.bwh.harvard.edu/cbase/CBaSE_v1.1.zip \
		-O ${tmpdir}/cbase.zip
	unzip -d ${tmpdir} ${tmpdir}/cbase.zip
	mv ${tmpdir}/CBaSE_v1.1/Auxiliary/*.gz ${DATASETS_CBASE}/
	mv ${tmpdir}/CBaSE_v1.1/Auxiliary/*.txt ${DATASETS_CBASE}/
	touch $@

ALL_DATASETS += $(CBASE_DATA)