#!/bin/bash

export INTOGEN_GENOME="hg19"
export INTOGEN_VEP="vep88"
export INTOGEN_RELEASE="develop"

nextflow run ../../intogen.nf -resume -profile local --input ./input --output ./output
