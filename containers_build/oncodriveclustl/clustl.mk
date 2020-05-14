# TODO use python as base image

CLUSTL_CONTAINER = $(INTOGEN_CONTAINERS)/oncodriveclustl.simg

clustl_container_srcdir = ${src_containers}/oncodriveclustl

clustl_container_src = ${clustl_container_srcdir}/Singularity

$(CLUSTL_CONTAINER): $(clustl_container_src) | $(INTOGEN_CONTAINERS)
	@echo Building OncodriveCLUSTL container
	cd ${clustl_container_srcdir} && \
		sudo singularity build ${tmpdir}/$(@F) Singularity
	mv ${tmpdir}/$(@F) $@
	sudo chown ${USER}: $@

CONTAINERS_SUDO += $(CLUSTL_CONTAINER)