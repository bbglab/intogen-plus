
FML_CONTAINER = $(INTOGEN_CONTAINERS)/oncodrivefml.simg

fml_container_srcdir = ${src_containers}/oncodrivefml

fml_container_src = ${fml_container_srcdir}/oncodrivefml_v2.conf \
	${fml_container_srcdir}/Singularity

$(FML_CONTAINER): $(fml_container_src) | $(INTOGEN_CONTAINERS)
	@echo Building OncodriveFML container
	${container_builder} ${fml_container_srcdir} $@

CONTAINERS_SUDO += $(FML_CONTAINER)