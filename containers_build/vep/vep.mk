
# Docs on the VEP docker image: https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html#docker

VEP_CONTAINER = $(INTOGEN_CONTAINERS)/vep.simg

vep_container_releases_file = ${src_containers}/vep/releases.txt
vep_container_version = `grep "^${ensembl}" ${vep_container_releases_file}`

$(VEP_CONTAINER): $(vep_container_releases_file) $$(ENSEMBL) | $(INTOGEN_CONTAINERS)
	@echo Building VEP container
	singularity build $@ docker://ensemblorg/ensembl-vep:release_${vep_container_version}

CONTAINERS_USER += $(VEP_CONTAINER)

# TODO remove
VEP_CONTAINER2 = $(INTOGEN_CONTAINERS)/vep${ensembl}.simg

vep_container_srcdir = ${src_containers}/vep

$(VEP_CONTAINER2): ${vep_container_srcdir}/Singularity.vep${ensembl} | $(INTOGEN_CONTAINERS)
	@echo Building VEP container 2
	cd ${vep_container_srcdir} && \
		sudo singularity build ${tmpdir}/$(@F) Singularity.vep${ensembl}
	mv ${tmpdir}/$(@F) $@
	sudo chown ${USER}: $@

CONTAINERS_SUDO += $(VEP_CONTAINER2)