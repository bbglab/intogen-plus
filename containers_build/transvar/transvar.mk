
IMAGE_TRANSVAR = $(CONTAINERS)/transvar.simg

$(IMAGE_TRANSVAR): | $(CONTAINERS)
	singularity build $@ docker://zhouwanding/transvar

ALL_TARGETS += $(IMAGE_TRANSVAR)