# TODO remove run script and use that in the nexflow file

SIGNATURE_CONTAINER = $(INTOGEN_CONTAINERS)/signature.simg

signature_container_srcdir = ${src_containers}/signature

signature_container_src = ${signature_container_srcdir}/run.sh \
	${signature_container_srcdir}/Singularity

$(SIGNATURE_CONTAINER): $(signature_container_src) | $(INTOGEN_CONTAINERS)
	@echo Building Signature container
	cd ${signature_container_srcdir} && \
		sudo singularity build ${tmpdir}/$(@F) Singularity
	mv ${tmpdir}/$(@F) $@
	sudo chown ${USER}: $@

CONTAINERS_SUDO += $(SIGNATURE_CONTAINER)