

CONTAINER_HOTMAPS = $(CONTAINERS)/hotmaps.simg

SRC_CONTAINER_HOTMAPS = ${CONTAINERS_SOURCE_FOLDER}/hotmaps

CONTAINER_HOTMAPS_FILES = $(wildcard ${SRC_CONTAINER_HOTMAPS}/*)

$(CONTAINER_HOTMAPS): ${CONTAINER_HOTMAPS_FILES} | $(CONTAINERS)
	@echo Building deconstructsig container
	cd ${SRC_CONTAINER_HOTMAPS} && \
		sudo singularity build $$(basename $@) Singularity
	mv ${SRC_CONTAINER_HOTMAPS}/$$(basename $@) $@
	sudo chown ${USER}: $@

TARGETS_CONTAINERS_SUDO += $(CONTAINER_HOTMAPS)