#!/bin/bash

export INTOGEN_GENOME="hg38"
export INTOGEN_VEP="vep92"
export INTOGEN_RELEASE="develop"

nextflow run ../../intogen.nf -resume -profile local --input ./input --output ./output
