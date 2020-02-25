
# TODO lowercase name
CONTAINER_MUTPANNING = $(CONTAINERS)/MutPanning.simg

SRC_CONTAINER_MUTPANNING = ${CONTAINERS_SOURCE_FOLDER}/mutpanning


$(CONTAINER_MUTPANNING): ${SRC_CONTAINER_MUTPANNING}/Singularity ${SRC_CONTAINER_MUTPANNING}/MutPanning.jar | $(CONTAINERS)
	@echo Building MutPanning container
	cd ${SRC_CONTAINER_MUTPANNING} && \
		sudo singularity build $$(basename $@) Singularity
	mv ${SRC_CONTAINER_MUTPANNING}/$$(basename $@) $@
	sudo chown ${USER}: $@

TARGETS_CONTAINERS_SUDO += $(CONTAINER_MUTPANNING)