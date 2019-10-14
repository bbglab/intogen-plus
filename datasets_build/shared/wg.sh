#!/bin/bash

set -e

if [ -z "${INTOGEN_RELEASE}" ]
then
      echo "ERROR: Define the INTOGEN_RELEASE variable"
      exit -1
fi


WG_REGIONS_FILE="../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/shared/wg.regions.gz"
python create_wg_regions.py ${INTOGEN_GENOME} 3 | gzip > ${WG_REGIONS_FILE}
