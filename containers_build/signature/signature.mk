
CONTAINER_SIGNATURE = $(CONTAINERS)/signature.simg

SRC_CONTAINER_SIGNATURE = ${CONTAINERS_SOURCE_FOLDER}/signature

$(CONTAINER_SIGNATURE): ${SRC_CONTAINER_SIGNATURE}/Singularity ${SRC_CONTAINER_SIGNATURE}/run.sh | $(CONTAINERS)
	@echo Building SIGNATURE container
	cd ${SRC_CONTAINER_SIGNATURE} && \
		sudo singularity build $$(basename $@) Singularity
	mv ${SRC_CONTAINER_SIGNATURE}/$$(basename $@) $@
	sudo chown ${USER}: $@

TARGETS_CONTAINERS_SUDO += $(CONTAINER_SIGNATURE)