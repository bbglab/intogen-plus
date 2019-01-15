#!/bin/bash

export INTOGEN_RELEASE="latest"
export INTOGEN_HOME=`echo "$(pwd)/../../" | xargs realpath`
nextflow run ${INTOGEN_HOME}/intogen.nf -resume -profile local --input ./input --output ./output