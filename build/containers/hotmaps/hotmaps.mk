
HOTMAPS_CONTAINER = $(INTOGEN_CONTAINERS)/hotmaps.simg

hotmaps_container_srcdir = ${src_containers}/hotmaps

hotmaps_container_src = $(wildcard ${hotmaps_container_srcdir}/*) \
	${hotmaps_container_srcdir}/hotmaps.sh ${hotmaps_container_srcdir}/Singularity

$(HOTMAPS_CONTAINER): $(hotmaps_container_src) | $(INTOGEN_CONTAINERS)
	@echo Building HotMAPS container
	${container_builder} ${hotmaps_container_srcdir} $@

CONTAINERS_SUDO += $(HOTMAPS_CONTAINER)