
SMREGIONS_CONTAINER = $(INTOGEN_CONTAINERS)/smregions.simg

smregions_container_srcdir = ${src_containers}/smregions

smregions_code_folder = ${src_containers}/../../smregions

smregions_container_src = $(wildcard ${smregions_code_folder}/*) \
	 $(wildcard ${smregions_code_folder}/smregions/*) \
	 ${smregions_container_srcdir}/smregions.conf \
	 ${smregions_container_srcdir}/Singularity

$(SMREGIONS_CONTAINER): $(smregions_container_src) | $(INTOGEN_CONTAINERS)
	@echo Building SMRegions container
	${container_builder} ${smregions_container_srcdir} $@

CONTAINERS_SUDO += $(SMREGIONS_CONTAINER)
