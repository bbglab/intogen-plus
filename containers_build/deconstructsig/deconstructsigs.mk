# TODO rename as deconstructsigs (last "s")

DECONSTRUCTSIGS_CONTAINER = $(INTOGEN_CONTAINERS)/deconstructsig.simg

deconstructsigs_container_srcdir = ${src_containers}/deconstructsig

deconstructsigs_container_src = $(wildcard ${deconstructsigs_container_srcdir}/*)

$(DECONSTRUCTSIGS_CONTAINER): $(deconstructsigs_container_src) | $(INTOGEN_CONTAINERS)
	@echo Building deconstructSigs container
	cd ${deconstructsigs_container_srcdir} && \
		sudo singularity build ${tmpdir}/$(@F) Singularity
	mv ${tmpdir}/$(@F) $@
	sudo chown ${USER}: $@

CONTAINERS_SUDO += $(DECONSTRUCTSIGS_CONTAINER)