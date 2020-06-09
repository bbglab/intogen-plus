
MUTPANNING_CONTAINER = $(INTOGEN_CONTAINERS)/mutpanning.simg

mutpanning_container_srcdir = ${src_containers}/mutpanning

mutpanning_container_src = ${mutpanning_container_srcdir}/MutPanning.jar \
	${mutpanning_container_srcdir}/Singularity

# FIXME remove jar file from repo and add instructions

$(MUTPANNING_CONTAINER): $(mutpanning_container_src) | $(INTOGEN_CONTAINERS)
	@echo Building MutPanning container
	${container_builder} ${mutpanning_container_srcdir} $@

CONTAINERS_SUDO += $(MUTPANNING_CONTAINER)