
CONTAINER_TRANSVAR = $(CONTAINERS)/transvar.simg

$(CONTAINER_TRANSVAR): | $(CONTAINERS)
	singularity build $@ docker://zhouwanding/transvar

TARGETS_CONTAINERS += $(CONTAINER_TRANSVAR)