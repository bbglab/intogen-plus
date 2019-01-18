#!/bin/bash

set -xe

# Configurable parameters
CADD_VERSION="v1.4"
CADD_GENOME=${INTOGEN_GENOME/hg/GRCh}
CADD_GENOME=${CADD_GENOME/GRCh19/GRCh37}

# We recomend to first download the whole CADD file and then use the file path
# instead of this URL (it will be faster and less likly to fail)
CADD_URL="https://krishna.gs.washington.edu/download/CADD/${CADD_VERSION}/${CADD_GENOME}/whole_genome_SNVs.tsv.gz"

# Build script
INTOGEN_DATASETS="../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}"
CDS_REGIONS="${INTOGEN_DATASETS}/shared/cds.regions.gz"
CADD_OUTPUT="${INTOGEN_DATASETS}/oncodrivefml/cadd.tsv.gz"

zcat ${CDS_REGIONS} | awk -v cadd="${CADD_URL}" '{system("tabix "cadd" "$1":"$2"-"$3)}' | gzip > ${CADD_OUTPUT}.tmp
zcat ${CADD_OUTPUT}.tmp | sort --parallel=8 -S 4G -k1,1 -k2,2n | uniq | bgzip > ${CADD_OUTPUT}
rm ${CADD_OUTPUT}.tmp
tabix -s 1 -b 2 -e 2 ${CADD_OUTPUT}
