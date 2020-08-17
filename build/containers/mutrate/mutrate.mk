
MUTRATE_CONTAINER = $(INTOGEN_CONTAINERS)/mutrate.simg

mutrate_container_srcdir = ${src_containers}/mutrate

mutrate_container_src = $(wildcard ${mutrate_container_srcdir}/*.py) \
	$(wildcard ${mutrate_container_srcdir}/*.sh) \
	${mutrate_container_srcdir}/Singularity


$(MUTRATE_CONTAINER): $(mutrate_container_src) | $(INTOGEN_CONTAINERS)
	@echo Building mutrate container
	${container_builder} ${hotmaps_container_srcdir} $@

CONTAINERS_SUDO += $(MUTRATE_CONTAINER)
