
DNDSCV_CONTAINER = $(INTOGEN_CONTAINERS)/dndscv.simg

dnds_container_srcdir = ${src_containers}/dndscv

dnds_container_src = $(wildcard ${dnds_container_srcdir}/*)

$(DNDSCV_CONTAINER): $(dnds_container_src) | $(INTOGEN_CONTAINERS)
	@echo Building dNdScv container
	cd ${dnds_container_srcdir} && \
		sudo singularity build ${tmpdir}/$(@F) Singularity
	mv ${tmpdir}/$(@F) $@
	sudo chown ${USER}: $@

CONTAINERS_SUDO += $(DNDSCV_CONTAINER)