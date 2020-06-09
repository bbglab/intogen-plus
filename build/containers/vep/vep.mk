
# Docs on the VEP docker image: https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html#docker

VEP_CONTAINER = $(INTOGEN_CONTAINERS)/vep.simg

vep_container_releases_file = ${src_containers}/vep/releases.txt
vep_container_version = `grep "^${ensembl}" ${vep_container_releases_file}`

$(VEP_CONTAINER): $(vep_container_releases_file) $$(ENSEMBL) | $(INTOGEN_CONTAINERS)
	@echo Building VEP container
	singularity build $@ docker://ensemblorg/ensembl-vep:release_${vep_container_version}

CONTAINERS_USER += $(VEP_CONTAINER)
