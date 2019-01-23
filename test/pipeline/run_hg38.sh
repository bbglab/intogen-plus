#!/bin/bash

export INTOGEN_GENOME="hg38"
export INTOGEN_VEP="vep92"
export INTOGEN_RELEASE="develop"
export INTOGEN_HOME=`echo "$(pwd)/../../" | xargs realpath`

nextflow run ${INTOGEN_HOME}/intogen.nf -resume -profile local --input ./input --output ./output
