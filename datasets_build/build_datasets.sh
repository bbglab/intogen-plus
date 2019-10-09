#!/bin/bash

set -e

export INTOGEN_GENOME=$1 "hg38"
export INTOGEN_VEP=$2 "vep92"
export INTOGEN_RELEASE=$3 "develop"

mkdir -p ../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}

# Shared
echo "Building shared datasets"
cd shared
./run.sh
cd ..

# Transvar
echo "Building Transvar datasets"
cd transvar
./build_hg38.sh
cd ..

# SmRegions
echo "Building SmRegions datasets"
cd smregions
./regions_pfam.sh
./panno.sh
cd ..

# dNdSCV
echo "Building dNdSCV datasets"
cd dndscv
./refcds_rda.sh
cd ..

# Mutrate
echo "Building Mutrate datasets"
cd mutrate
python cosmic_to_exome.py
cd ..

# OncodriveFML
echo "Building OncodriveFML datasets"
cd oncodrivefml
./cadd.sh
cd ..

# PTMS
echo "Building PTMA datasets"
cd ptms
./run_ptms.sh
cd ..

# boostDM
echo "Building boostDM datasets"
cd boostDM
./boostDM_features.sh
cd ..

# Combination
echo "Building combination datasets"
cd combination
./run_negative_set.sh
./run_cgc.sh
cd ..

sudo chown -R $USER: ../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}

