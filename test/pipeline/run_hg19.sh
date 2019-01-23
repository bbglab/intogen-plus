#!/bin/bash

export INTOGEN_GENOME="hg19"
export INTOGEN_VEP="vep88"
export INTOGEN_RELEASE="develop"
export INTOGEN_HOME=`echo "$(pwd)/../../" | xargs realpath`

nextflow run ${INTOGEN_HOME}/intogen.nf -resume -profile local --input ./input --output ./output
