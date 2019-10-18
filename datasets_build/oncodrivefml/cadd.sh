#!/bin/bash

set -xe

if [ -z "${INTOGEN_RELEASE}" ]
then
      echo "ERROR: Define the INTOGEN_RELEASE variable"
      exit -1
fi

if [ -z "${INTOGEN_VEP}" ]
then
      echo "ERROR: Define the INTOGEN_VEP variable"
      exit -1
fi

if [ -z "${INTOGEN_GENOME}" ]
then
      echo "ERROR: Define the INTOGEN_GENOME variable"
      exit -1
fi

# Configurable parameters
CADD_VERSION="v1.4"
CADD_GENOME=${INTOGEN_GENOME/hg/GRCh}
CADD_GENOME=${CADD_GENOME/GRCh19/GRCh37}

# We recomend to first download the whole CADD file with `wget -c ${CADD_URL}`and then
# use the file path instead of this URL (it will be faster and less likely to fail).
# Also note that in order to use tabix with a remote download, HTTP should be used intead of HTTPS
CADD_URL="http://krishna.gs.washington.edu/download/CADD/${CADD_VERSION}/${CADD_GENOME}/whole_genome_SNVs.tsv.gz"

# Build script
INTOGEN_DATASETS="../../datasets/"${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}
mkdir -p "${INTOGEN_DATASETS}/oncodrivefml"
CDS_REGIONS="${INTOGEN_DATASETS}/shared/cds.regions.gz"
CADD_OUTPUT="${INTOGEN_DATASETS}/oncodrivefml/cadd.tsv.gz"

zcat ${CDS_REGIONS} | tail -n +2| awk -v cadd="${CADD_URL}" '{system("tabix "cadd" "$1":"$2"-"$3)}' | gzip > ${CADD_OUTPUT}.tmp
zcat ${CADD_OUTPUT}.tmp | sort --parallel=8 -S 4G -k1,1 -k2,2n | uniq | bgzip > ${CADD_OUTPUT}
rm ${CADD_OUTPUT}.tmp
tabix -s 1 -b 2 -e 2 ${CADD_OUTPUT}
