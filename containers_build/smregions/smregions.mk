# TODO remove cpus from config file, or maybe move config file to datasets

SMREGIONS_CONTAINER = $(CONTAINERS)/smregions.simg

smregions_container_srcdir = ${src_containers}/smregions

smregions_container_src = $(wildcard ${smregions_container_srcdir}/src/*) \
	 $(wildcard ${smregions_container_srcdir}/src/smregions/*) \
	 ${smregions_container_srcdir}/smregions.conf \
	 ${smregions_container_srcdir}/Singularity

$(SMREGIONS_CONTAINER): $(smregions_container_src) | $(CONTAINERS)
	@echo Building SMRegions container
	cd ${smregions_container_srcdir} && \
		sudo singularity build ${tmpdir}/$(@F) Singularity
	mv ${tmpdir}/$(@F) $@
	sudo chown ${USER}: $@

CONTAINERS_SUDO += $(SMREGIONS_CONTAINER)