
#CONTAINER_TRANSVAR = $(CONTAINERS)/transvar.simg
#
#$(CONTAINER_TRANSVAR): | $(CONTAINERS)
#	singularity build $@ docker://zhouwanding/transvar
#
#TARGETS_CONTAINERS += $(CONTAINER_TRANSVAR)


CONTAINER_TRANSVAR = $(CONTAINERS)/transvar.simg

SRC_CONTAINER_TRANSVAR = ${CONTAINERS_SOURCE_FOLDER}/transvar

$(CONTAINER_TRANSVAR): ${SRC_CONTAINER_TRANSVAR}/Singularity | $(CONTAINERS)
	@echo Building TRANSVAR container
	cd ${SRC_CONTAINER_TRANSVAR} && \
		sudo singularity build $$(basename $@) Singularity
	mv ${SRC_CONTAINER_TRANSVAR}/$$(basename $@) $@
	sudo chown ${USER}: $@

TARGETS_CONTAINERS_SUDO += $(CONTAINER_TRANSVAR)