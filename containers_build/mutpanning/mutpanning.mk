
MUTPANNING_CONTAINER = $(INTOGEN_CONTAINERS)/mutpanning.simg

mutpanning_container_srcdir = ${src_containers}/mutpanning

mutpanning_container_src = ${mutpanning_container_srcdir}/MutPanning.jar \
	${mutpanning_container_srcdir}/Singularity

# TODO remove jar file from repo and add instructions

$(MUTPANNING_CONTAINER): $(mutpanning_container_src) | $(INTOGEN_CONTAINERS)
	@echo Building MutPanning container
	cd ${mutpanning_container_srcdir} && \
		sudo singularity build ${tmpdir}/$(@F) Singularity
	mv ${tmpdir}/$(@F) $@
	sudo chown ${USER}: $@

CONTAINERS_SUDO += $(MUTPANNING_CONTAINER)