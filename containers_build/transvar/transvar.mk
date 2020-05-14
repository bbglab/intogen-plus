
#$(CONTAINER_TRANSVAR): | $(INTOGEN_CONTAINERS)
#	singularity build $@ docker://zhouwanding/transvar


TRANSVAR_CONTAINER = $(INTOGEN_CONTAINERS)/transvar.simg

transvar_container_srcdir = ${src_containers}/transvar

transvar_container_src = ${transvar_container_srcdir}/Singularity

$(TRANSVAR_CONTAINER): $(transvar_container_src) | $(INTOGEN_CONTAINERS)
	@echo Building TransVar container
	cd ${transvar_container_srcdir} && \
		sudo singularity build ${tmpdir}/$(@F) Singularity
	mv ${tmpdir}/$(@F) $@
	sudo chown ${USER}: $@

CONTAINERS_SUDO += $(TRANSVAR_CONTAINER)