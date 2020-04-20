
cbase_dir = $(DATASETS)/cbase

CBASE_DATA = ${cbase_dir}/.checkpoint

$(CBASE_DATA):
	@echo Building CBaSE datasets
	mkdir -p ${cbase_dir}
	wget -c http://genetics.bwh.harvard.edu/cbase/CBaSE_v1.1.zip \
		-O ${tmpdir}/cbase.zip
	unzip -d ${tmpdir} ${tmpdir}/cbase.zip
	mv ${tmpdir}/CBaSE_v1.1/Auxiliary/*.gz ${cbase_dir}/
	mv ${tmpdir}/CBaSE_v1.1/Auxiliary/*.txt ${cbase_dir}/
	touch $@


DATASETS_TARGETS += $(CBASE_DATA)