

CONTAINER_MUTRATE = $(CONTAINERS)/mutrate.simg

SRC_CONTAINER_MUTRATE = ${CONTAINERS_SOURCE_FOLDER}/mutrate

CONTAINER_MUTRATE_FILES = $(wildcard ${SRC_CONTAINER_MUTRATE}/*.py) $(wildcard ${SRC_CONTAINER_MUTRATE}/*.sh)


$(CONTAINER_MUTRATE): ${SRC_CONTAINER_MUTRATE}/Singularity ${CONTAINER_MUTRATE_FILES} | $(CONTAINERS)
	@echo Building mutrate container
	cd ${SRC_CONTAINER_MUTRATE} && \
		sudo singularity build $$(basename $@) Singularity
	mv ${SRC_CONTAINER_MUTRATE}/$$(basename $@) $@
	sudo chown ${USER}: $@

TARGETS_CONTAINERS_SUDO += $(CONTAINER_MUTRATE)