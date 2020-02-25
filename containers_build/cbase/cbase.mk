

CONTAINER_CBASE = $(CONTAINERS)/cbase.simg

SRC_CONTAINER_CBASE = ${CONTAINERS_SOURCE_FOLDER}/cbase

$(CONTAINER_CBASE): ${SRC_CONTAINER_CBASE}/Singularity ${SRC_CONTAINER_CBASE}/cbase.py | $(CONTAINERS)
	@echo Building CBaSE container
	cd ${SRC_CONTAINER_CBASE} && \
		sudo singularity build $$(basename $@) Singularity
	mv ${SRC_CONTAINER_CBASE}/$$(basename $@) $@
	sudo chown ${USER}: $@

TARGETS_CONTAINERS_SUDO += $(CONTAINER_CBASE)