
CLUSTL_CONTAINER = $(INTOGEN_CONTAINERS)/oncodriveclustl.simg

clustl_container_srcdir = ${src_containers}/oncodriveclustl

clustl_container_src = ${clustl_container_srcdir}/Singularity

$(CLUSTL_CONTAINER): $(clustl_container_src) | $(INTOGEN_CONTAINERS)
	@echo Building OncodriveCLUSTL container
	${container_builder} ${clustl_container_srcdir} $@

CONTAINERS_SUDO += $(CLUSTL_CONTAINER)