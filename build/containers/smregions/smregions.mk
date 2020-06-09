
SMREGIONS_CONTAINER = $(INTOGEN_CONTAINERS)/smregions.simg

smregions_container_srcdir = ${src_containers}/smregions

smregions_container_src = ${smregions_container_srcdir}/smregions.conf \
	 ${smregions_container_srcdir}/Singularity

$(SMREGIONS_CONTAINER): $(smregions_container_src) | $(INTOGEN_CONTAINERS)
	@echo Building SMRegions container
	${container_builder} ${smregions_container_srcdir} $@


CONTAINERS_SUDO += $(SMREGIONS_CONTAINER)