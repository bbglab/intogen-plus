
FML_CONTAINER = $(INTOGEN_CONTAINERS)/oncodrivefml.simg

fml_container_srcdir = ${src_containers}/oncodrivefml

# TODO put config file as dataset
# TODO make reference genome as variable in the config

fml_container_src = ${fml_container_srcdir}/oncodrivefml_v2.conf \
	${fml_container_srcdir}/Singularity

$(FML_CONTAINER): $(fml_container_src) | $(INTOGEN_CONTAINERS)
	@echo Building OncodriveFML container
	cd ${fml_container_srcdir} && \
		sudo singularity build ${tmpdir}/$(@F) Singularity
	mv ${tmpdir}/$(@F) $@
	sudo chown ${USER}: $@

CONTAINERS_SUDO += $(FML_CONTAINER)