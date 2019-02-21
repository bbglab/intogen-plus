#!/usr/bin/env bash

set -e

if [ -z "${INTOGEN_RELEASE}" ]
then
      echo "ERROR: Define the INTOGEN_RELEASE variable"
      exit -1
fi


CDS_REGIONS_FILE="../../datasets/${INTOGEN_RELEASE}/shared/cds.regions.gz"
TRIPLETS_FILE="../../datasets/${INTOGEN_RELEASE}/shared/triplets.gz"
CONSEQUENCES_FILE="../../datasets/${INTOGEN_RELEASE}/shared/consequences.gz"

python count.py -r ${CDS_REGIONS_FILE} -t ${TRIPLETS_FILE} -c ${CONSEQUENCES_FILE} -v 92 -g hg38
