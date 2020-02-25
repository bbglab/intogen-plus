
IMAGE_VEP = $(CONTAINERS)/vep.simg

VEP_RELEASES_FILE = ${CONTAINERS_SOURCE_FOLDER}/vep/releases.txt
VEP_RELEASE_VERSION = `grep "^${ENSEMBL}" ${VEP_RELEASES_FILE}`

$(IMAGE_VEP): ${VEP_RELEASES_FILE} | $(CONTAINERS)
	singularity build $@ docker://ensemblorg/ensembl-vep:release_${VEP_RELEASE_VERSION}

ALL_TARGETS += $(IMAGE_VEP)