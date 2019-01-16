#!/bin/bash

set -e

export INTOGEN_VEP="vep92"
export INTOGEN_GENOME="hg38"
export INTOGEN_RELEASE="develop"

mkdir -p ../containers/${INTOGEN_RELEASE}

# CBase
echo "Building cbase.simg"
cd cbase
sudo singularity build ../../containers/${INTOGEN_RELEASE}/cbase.simg Singularity
cd ..

# Combination
echo "Building combination.simg"
cd combination
sudo singularity build ../../containers/${INTOGEN_RELEASE}/combination.simg Singularity
cd ..

# DeconstructSig
echo "Building deconstructsig.simg"
cd deconstructsig
sudo singularity build ../../containers/${INTOGEN_RELEASE}/deconstructsig.simg Singularity
cd ..

# DndsCV
echo "Building dndscv.simg"
cd dndscv
sudo singularity build ../../containers/${INTOGEN_RELEASE}/dndscv.simg Singularity
cd ..

# Hotmaps
echo "Building hotmaps.simg"
cd hotmaps
sudo singularity build ../../containers/${INTOGEN_RELEASE}/hotmaps.simg Singularity
cd ..

# Mutrate
echo "Building mutrate.simg"
cd mutrate
sudo singularity build ../../containers/${INTOGEN_RELEASE}/mutrate.simg Singularity
cd ..

# OncodriveCLUSTL
echo "Building oncodriveclustl.simg"
cd oncodriveclustl
sudo singularity build ../../containers/${INTOGEN_RELEASE}/oncodriveclustl.simg Singularity
cd ..

# OncodriveFML
echo "Building oncodrivefml.simg"
cd oncodrivefml
sudo singularity build ../../containers/${INTOGEN_RELEASE}/oncodrivefml.simg Singularity
cd ..

# SmRegions 
echo "Building smregions.simg"
cd smregions
sudo singularity build ../../containers/${INTOGEN_RELEASE}/smregions.simg Singularity
cd ..

# Variant Effect Predictor 
echo "Building ${INTOGEN_VEP}.simg"
cd vep
sudo singularity build ../../containers/${INTOGEN_RELEASE}/${INTOGEN_VEP}.simg Singularity.${INTOGEN_VEP}
cd ..
