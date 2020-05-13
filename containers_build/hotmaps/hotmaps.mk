
HOTMAPS_CONTAINER = $(CONTAINERS)/hotmaps.simg

hotmaps_container_srcdir = ${src_containers}/hotmaps

hotmaps_container_src = $(wildcard ${hotmaps_container_srcdir}/*)

$(HOTMAPS_CONTAINER): $(hotmaps_container_src) | $(CONTAINERS)
	@echo Building HotMAPS container
	cd ${hotmaps_container_srcdir} && \
		sudo singularity build ${tmpdir}/$(@F) Singularity
	mv ${tmpdir}/$(@F) $@
	sudo chown ${USER}: $@

CONTAINERS_SUDO += $(HOTMAPS_CONTAINER)