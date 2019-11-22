#!/usr/bin/env bash

set -e

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


CDS_REGIONS_FILE="../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/shared/cds.regions.gz"
WG_REGIONS_FILE="../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/shared/wg.regions.gz"
COUNT_CDS="../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/shared/cds.counts.gz"
COUNT_WG="../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/shared/wg.counts.gz"
CONSEQUENCE_CDS="../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/shared/consequences.pickle.gz"
TRIPLETS_CDS="../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/shared/triplets.json.gz"

bgsignature count -r ${CDS_REGIONS_FILE} -s 3 -g ${INTOGEN_GENOME} --cores 12 --collapse --exclude-N -o ${COUNT_CDS}
bgsignature count -r ${WG_REGIONS_FILE} -s 3 -g ${INTOGEN_GENOME} --cores 12 --collapse --exclude-N -o ${COUNT_WG}
python count.py -r ${CDS_REGIONS_FILE} -t ${TRIPLETS_CDS} -c ${CONSEQUENCE_CDS} -v ${INTOGEN_VEP#"vep"} -g ${INTOGEN_GENOME}