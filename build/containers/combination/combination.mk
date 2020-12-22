
COMBINATION_CONTAINER = $(INTOGEN_CONTAINERS)/intogen-combination.simg

combination_container_srcdir = ${src_containers}/combination

combination_code_folder = ${src_containers}/../../combination

combination_container_src = ${combination_container_srcdir}/Singularity \
	${combination_code_folder}/setup.py \
	$(wildcard ${combination_code_folder}/intogen_combination/*) \
	$(wildcard ${combination_code_folder}/*)

$(COMBINATION_CONTAINER): ${combination_container_src} | $(INTOGEN_CONTAINERS)
	@echo Building combination container
	${container_builder} ${combination_container_srcdir} $@

CONTAINERS_SUDO += $(COMBINATION_CONTAINER)