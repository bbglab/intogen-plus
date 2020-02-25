
CONTAINER_CLUSTL = $(CONTAINERS)/oncodriveclustl.simg

SRC_CONTAINER_CLUSTL = ${CONTAINERS_SOURCE_FOLDER}/oncodriveclustl

$(CONTAINER_CLUSTL): ${SRC_CONTAINER_CLUSTL}/Singularity | $(CONTAINERS)
	@echo Building CLUSTL container
	cd ${SRC_CONTAINER_CLUSTL} && \
		sudo singularity build $$(basename $@) Singularity
	mv ${SRC_CONTAINER_CLUSTL}/$$(basename $@) $@
	sudo chown ${USER}: $@

TARGETS_CONTAINERS_SUDO += $(CONTAINER_CLUSTL)