
CORE_CONTAINER = $(INTOGEN_CONTAINERS)/intogen-core.simg

core_container_srcdir = ${src_containers}/core

core_code_folder = ${src_containers}/../../core

core_container_src = ${core_container_srcdir}/Singularity \
	${core_container_srcdir}/get_field.sh \
	${core_code_folder}/setup.py \
	$(wildcard ${core_code_folder}/intogen_core/*.py) \
	$(wildcard ${core_code_folder}/intogen_core/*/*.py) \
	$(wildcard ${core_code_folder}/intogen_core/*/*/*.py)

$(CORE_CONTAINER): ${core_container_src} | $(INTOGEN_CONTAINERS)
	@echo Building core container
	${container_builder} ${core_container_srcdir} $@

CONTAINERS_SUDO += $(CORE_CONTAINER)