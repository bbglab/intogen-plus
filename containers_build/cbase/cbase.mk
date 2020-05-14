

CBASE_CONTAINER = $(INTOGEN_CONTAINERS)/cbase.simg

cbase_container_srcdir = ${src_containers}/cbase

cbase_container_src = ${cbase_container_srcdir}/cbase.py \
	${cbase_container_srcdir}/Singularity

$(CBASE_CONTAINER): $(cbase_container_src) | $(INTOGEN_CONTAINERS)
	@echo Building CBaSE container
	cd ${cbase_container_srcdir} && \
		sudo singularity build ${tmpdir}/$(@F) Singularity
	mv ${tmpdir}/$(@F) $@
	sudo chown ${USER}: $@

CONTAINERS_SUDO += $(CBASE_CONTAINER)