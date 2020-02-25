

CONTAINER_COMBINATION = $(CONTAINERS)/combination.simg

SRC_CONTAINER_COMBINATION = ${CONTAINERS_SOURCE_FOLDER}/combination

CONTAINER_COMBINATION_FILES = $(wildcard ${SRC_CONTAINER_COMBINATION}/src/*) \
	 $(wildcard ${SRC_CONTAINER_COMBINATION}/src/config/*) \
	 $(wildcard ${SRC_CONTAINER_COMBINATION}/src/evaluation/*) \
	 $(wildcard ${SRC_CONTAINER_COMBINATION}/src/qc/*)


$(CONTAINER_COMBINATION): ${SRC_CONTAINER_COMBINATION}/Singularity ${CONTAINER_COMBINATION_FILES} | $(CONTAINERS)
	@echo Building combination container
	cd ${SRC_CONTAINER_COMBINATION} && \
		sudo singularity build $$(basename $@) Singularity
	mv ${SRC_CONTAINER_COMBINATION}/$$(basename $@) $@
	sudo chown ${USER}: $@

TARGETS_CONTAINERS_SUDO += $(CONTAINER_COMBINATION)