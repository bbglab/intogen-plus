

CONTAINER_DNDSCV = $(CONTAINERS)/dndscv.simg

SRC_CONTAINER_DNDSCV = ${CONTAINERS_SOURCE_FOLDER}/dndscv

CONTAINER_DNDSCV_FILES = $(wildcard ${SRC_CONTAINER_DNDSCV}/*)

$(CONTAINER_DNDSCV): ${CONTAINER_DNDSCV_FILES} | $(CONTAINERS)
	@echo Building deconstructsig container
	cd ${SRC_CONTAINER_DNDSCV} && \
		sudo singularity build $$(basename $@) Singularity
	mv ${SRC_CONTAINER_DNDSCV}/$$(basename $@) $@
	sudo chown ${USER}: $@

TARGETS_CONTAINERS_SUDO += $(CONTAINER_DNDSCV)