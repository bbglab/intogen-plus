

COMBINATION_CONTAINER = $(CONTAINERS)/combination.simg

combination_container_srcdir = ${src_containers}/combination

combination_container_src = $(wildcard ${combination_container_srcdir}/src/*) \
	 $(wildcard ${combination_container_srcdir}/src/config/*) \
	 $(wildcard ${combination_container_srcdir}/src/evaluation/*) \
	 $(wildcard ${combination_container_srcdir}/src/qc/*) \
	 ${combination_container_srcdir}/Singularity


$(COMBINATION_CONTAINER): $(combination_container_src) | $(CONTAINERS)
	@echo Building combination container
	cd ${combination_container_srcdir} && \
		sudo singularity build ${tmpdir}/$(@F) Singularity
	mv ${tmpdir}/$(@F) $@
	sudo chown ${USER}: $@

CONTAINERS_SUDO += $(COMBINATION_CONTAINER)