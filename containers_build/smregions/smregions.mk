

CONTAINER_SMREGIONS = $(CONTAINERS)/smregions.simg

SRC_CONTAINER_SMREGIONS = ${CONTAINERS_SOURCE_FOLDER}/smregions

CONTAINER_SMREGIONS_FILES = $(wildcard ${SRC_CONTAINER_SMREGIONS}/src/*) \
	 $(wildcard ${SRC_CONTAINER_SMREGIONS}/src/smregions/*) \
	 ${SRC_CONTAINER_SMREGIONS}/smregions.conf

$(CONTAINER_SMREGIONS): ${SRC_CONTAINER_SMREGIONS}/Singularity ${CONTAINER_SMREGIONS_FILES} | $(CONTAINERS)
	@echo Building SMREGIONS container
	cd ${SRC_CONTAINER_SMREGIONS} && \
		sudo singularity build $$(basename $@) Singularity
	mv ${SRC_CONTAINER_SMREGIONS}/$$(basename $@) $@
	sudo chown ${USER}: $@

TARGETS_CONTAINERS_SUDO += $(CONTAINER_SMREGIONS)