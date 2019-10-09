#!/bin/bash

set -e

export INTOGEN_GENOME=$1 "hg38"
export INTOGEN_VEP=$2 "vep92"
export INTOGEN_RELEASE=$3  "v20191009"
export INTOGEN_CONTAINERS=$4 #"/workspace/projects/intogen_2017/pipeline/containers/"
mkdir -p ../containers/${INTOGEN_RELEASE}

# CBase

if [ ! -f ${INTOGEN_CONTAINERS}/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/cbase.simg ]; then
    echo "Building cbase.simg..."
    cd cbase
    sudo singularity build ../../containers/${INTOGEN_RELEASE}/cbase.simg Singularity
    mv ../../containers/${INTOGEN_RELEASE}/cbase.simg ${INTOGEN_CONTAINERS}/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/
    cd ..
fi


# Combination
if [ ! -f ${INTOGEN_CONTAINERS}/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/combination.simg ]; then
    echo "Building combination.simg"
    cd combination
    sudo singularity build ../../containers/${INTOGEN_RELEASE}/combination.simg Singularity
    mv ../../containers/${INTOGEN_RELEASE}/combination.simg ${INTOGEN_CONTAINERS}/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/
    cd ..
fi


# DeconstructSig

if [ ! -f ${INTOGEN_CONTAINERS}/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/deconstructsig.simg ]; then
    echo "Building deconstructsig.simg"
    cd deconstructsig
    sudo singularity build ../../containers/${INTOGEN_RELEASE}/deconstructsig.simg Singularity
    mv ../../containers/${INTOGEN_RELEASE}/deconstructsig.simg ${INTOGEN_CONTAINERS}/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/
    cd ..
fi


# DndsCV
if [ ! -f ${INTOGEN_CONTAINERS}/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/dndscv.simg ]; then
    echo "Building dndscv.simg"
    cd dndscv
    sudo singularity build ../../containers/${INTOGEN_RELEASE}/dndscv.simg Singularity
    mv ../../containers/${INTOGEN_RELEASE}/dndscv.simg ${INTOGEN_CONTAINERS}/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/
    cd ..
fi

# Hotmaps
if [ ! -f ${INTOGEN_CONTAINERS}/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/hotmaps.simg ]; then
    echo "Building hotmaps.simg"
    cd hotmaps
    sudo singularity build ../../containers/${INTOGEN_RELEASE}/hotmaps.simg Singularity
    mv ../../containers/${INTOGEN_RELEASE}/hotmaps.simg ${INTOGEN_CONTAINERS}/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/
    cd ..
fi

# Mutrate
if [ ! -f ${INTOGEN_CONTAINERS}/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/mutrate.simg ]; then
    echo "Building mutrate.simg"
    cd mutrate
    sudo singularity build ../../containers/${INTOGEN_RELEASE}/mutrate.simg Singularity
    mv ../../containers/${INTOGEN_RELEASE}/mutrate.simg ${INTOGEN_CONTAINERS}/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/
    cd ..
fi

# OncodriveCLUSTL
if [ ! -f ${INTOGEN_CONTAINERS}/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/oncodriveclustl.simg ]; then
    echo "Building oncodriveclustl.simg"
    cd oncodriveclustl
    sudo singularity build ../../containers/${INTOGEN_RELEASE}/oncodriveclustl.simg Singularity
    mv ../../containers/${INTOGEN_RELEASE}/oncodriveclustl.simg ${INTOGEN_CONTAINERS}/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/
    cd ..
fi


# OncodriveFML
if [ ! -f ${INTOGEN_CONTAINERS}/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/oncodrivefml.simg ]; then
    echo "Building oncodrivefml.simg"
    cd oncodrivefml
    sudo singularity build ../../containers/${INTOGEN_RELEASE}/oncodrivefml.simg Singularity
    mv ../../containers/${INTOGEN_RELEASE}/oncodrivefml.simg ${INTOGEN_CONTAINERS}/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/
    cd ..
fi


# SmRegions
if [ ! -f ${INTOGEN_CONTAINERS}/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/smregions.simg ]; then
    echo "Building smregions.simg"
    cd smregions
    sudo singularity build ../../containers/${INTOGEN_RELEASE}/smregions.simg Singularity
    mv ../../containers/${INTOGEN_RELEASE}/smregions.simg ${INTOGEN_CONTAINERS}/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/
    cd ..
fi

# Variant Effect Predictor
if [ ! -f ${INTOGEN_CONTAINERS}/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/${INTOGEN_VEP}.simg ]; then
    echo "Building ${INTOGEN_VEP}.simg"
    cd vep
    sudo singularity build ../../containers/${INTOGEN_RELEASE}/${INTOGEN_VEP}.simg Singularity.${INTOGEN_VEP}
    mv ../../containers/${INTOGEN_RELEASE}/${INTOGEN_VEP}.simg ${INTOGEN_CONTAINERS}/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/
    cd ..
fi

# Transvar
if [ ! -f ${INTOGEN_CONTAINERS}/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/transvar.simg ]; then
    echo "Building transvar.simg"
    cd transvar
    sudo singularity build ../../containers/${INTOGEN_RELEASE}/transvar.simg Singularity
    mv ../../containers/${INTOGEN_RELEASE}/transvar.simg ${INTOGEN_CONTAINERS}/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/
    cd ..
fi
sudo chown -R $USER: ../containers/${INTOGEN_RELEASE}/*.simg

