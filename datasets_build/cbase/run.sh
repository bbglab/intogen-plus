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


wget -c http://genetics.bwh.harvard.edu/cbase/CBaSE_v1.1.zip
unzip CBaSE_v1.1.zip
cp CBaSE_v1.1/Auxiliary/*.gz ../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/cbase/.
cp CBaSE_v1.1/Auxiliary/*.txt ../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/cbase/.
rm -r CBaSE_v1.1/