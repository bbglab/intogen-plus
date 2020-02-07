
FOLDER ?= cbase

.PHONY: all clean

CHECKPOINT = ${FOLDER}/.checkpoint

all: $(CHECKPOINT)

$(CHECKPOINT):
	@echo Building CBaSE datasets
	mkdir -p ${FOLDER}
	wget -c http://genetics.bwh.harvard.edu/cbase/CBaSE_v1.1.zip \
		-O ${FOLDER}/cbase.zip
	unzip -d ${FOLDER} ${FOLDER}/cbase.zip
	mv ${FOLDER}/CBaSE_v1.1/Auxiliary/*.gz ${FOLDER}/
	mv ${FOLDER}/CBaSE_v1.1/Auxiliary/*.txt ${FOLDER}/
	rm -r ${FOLDER}/__MACOSX/
	rm ${FOLDER}/cbase.zip
	touch $@


clean:
	rm -f $(CHECKPOINT)