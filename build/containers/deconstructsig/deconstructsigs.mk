
DECONSTRUCTSIGS_CONTAINER = $(INTOGEN_CONTAINERS)/deconstructsigs.simg

deconstructsigs_container_srcdir = ${src_containers}/deconstructsig

deconstructsigs_container_src = $(wildcard ${deconstructsigs_container_srcdir}/*)

$(DECONSTRUCTSIGS_CONTAINER): $(deconstructsigs_container_src) | $(INTOGEN_CONTAINERS)
	@echo Building deconstructSigs container
	${container_builder} ${deconstructsigs_container_srcdir} $@

CONTAINERS_SUDO += $(DECONSTRUCTSIGS_CONTAINER)