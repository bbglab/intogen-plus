
# Docs on the VEP docker image: https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html#docker

CONTAINER_VEP = $(CONTAINERS)/vep.simg

VEP_RELEASES_FILE = ${CONTAINERS_SOURCE_FOLDER}/vep/releases.txt
VEP_RELEASE_VERSION = `grep "^${ENSEMBL}" ${VEP_RELEASES_FILE}`

$(CONTAINER_VEP): ${VEP_RELEASES_FILE} | $(CONTAINERS)
	singularity build $@ docker://ensemblorg/ensembl-vep:release_${VEP_RELEASE_VERSION}

TARGETS_CONTAINERS += $(CONTAINER_VEP)