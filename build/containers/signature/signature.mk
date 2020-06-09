
SIGNATURE_CONTAINER = $(INTOGEN_CONTAINERS)/signature.simg

signature_container_srcdir = ${src_containers}/signature

signature_container_src = ${signature_container_srcdir}/Singularity

$(SIGNATURE_CONTAINER): ${signature_container_src} | $(INTOGEN_CONTAINERS)
	@echo Building bgSignature container
	${container_builder} ${signature_container_srcdir} $@

CONTAINERS_SUDO += $(SIGNATURE_CONTAINER)