
seismic_CONTAINER = $(INTOGEN_CONTAINERS)/seismic.simg

seismic_container_srcdir = ${src_containers}/seismic

seismic_container_src = ${seismic_container_srcdir}/Singularity

$(seismic_CONTAINER): ${seismic_container_src} | $(INTOGEN_CONTAINERS)
	@echo Building seismic container
	${container_builder} ${seismic_container_srcdir} $@

CONTAINERS_SUDO += $(seismic_CONTAINER)