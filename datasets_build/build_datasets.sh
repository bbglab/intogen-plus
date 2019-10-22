#!/bin/bash

set -ex

export INTOGEN_GENOME=$1 "hg38"
export INTOGEN_VEP=$2 "vep92"
export INTOGEN_RELEASE=$3 "v20191009"

mkdir -p ../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}

# Bgdata
echo "- Bgdata datasets"
BGDATA_LOCAL="../datasets/"${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}"/bgdata"
BGDATA_OFFLINE="FALSE"
export BGDATA_LOCAL=${BGDATA_LOCAL}
export BGDATA_OFFLINE=${BGDATA_OFFLINE}
bgdata get datasets/genomereference/${INTOGEN_GENOME}
# Redo the download (it will not download the data again) to force to get the hg38.master file
bgdata get datasets/genomereference/${INTOGEN_GENOME}

# Shared
echo "- Building shared datasets"
cd shared
./run.sh
cd ..

# Transvar
echo "- Building Transvar datasets"
cd transvar
./build_hg38.sh
cd ..

# SmRegions
echo "- Building SmRegions datasets"
cd smregions
./regions_pfam.sh
cd ..

# dNdSCV
echo "- Building dNdSCV datasets"
cd dndscv
./refcds_rda.sh
cd ..

# Mutrate
echo "- Building Mutrate datasets"    # FIXME: failing
cd mutrate
python cosmic_to_exome.py
cd ..

# OncodriveFML
echo "- Building OncodriveFML datasets"
cd oncodrivefml
./cadd.sh
cd ..

# PTMS
echo "- Building PTMA datasets"
#cd ptms
#./run_ptms.sh
#cd ..

# boostDM
echo "- Building boostDM datasets"
#cd boostDM
#./boostDM_features.sh
#cd ..

# Combination
echo "- Building combination datasets"
cd combination
./run.sh
cd ..

chown -R $USER: ../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}
