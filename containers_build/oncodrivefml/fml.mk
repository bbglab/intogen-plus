
CONTAINER_FML = $(CONTAINERS)/oncodrivefml.simg

SRC_CONTAINER_FML = ${CONTAINERS_SOURCE_FOLDER}/oncodrivefml

$(CONTAINER_FML): ${SRC_CONTAINER_FML}/Singularity ${SRC_CONTAINER_FML}/oncodrivefml_v2.conf | $(CONTAINERS)
	@echo Building FML container
	cd ${SRC_CONTAINER_FML} && \
		sudo singularity build $$(basename $@) Singularity
	mv ${SRC_CONTAINER_FML}/$$(basename $@) $@
	sudo chown ${USER}: $@

TARGETS_CONTAINERS_SUDO += $(CONTAINER_FML)