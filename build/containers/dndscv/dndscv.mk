
DNDSCV_CONTAINER = $(INTOGEN_CONTAINERS)/dndscv.simg

dndscv_container_srcdir = ${src_containers}/dndscv

dndscv_container_src = ${dndscv_container_srcdir}/dndscv.R \
	${dndscv_container_srcdir}/dndscv.tar.gz \
	${dndscv_container_srcdir}/Singularity

$(DNDSCV_CONTAINER): $(dndscv_container_src) | $(INTOGEN_CONTAINERS)
	@echo Building dNdScv container
	${container_builder} ${dndscv_container_srcdir} $@

CONTAINERS_SUDO += $(DNDSCV_CONTAINER)