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

echo -e "\tBuilding ensembl canonical transcripts"
./ensembl_canonical_transcripts.sh

echo -e "\tBuilding CDS annotations"
./cds.sh

echo -e "\tBuilding WG annotations"
./wg.sh

echo -e "\tBuilding counts"
./counts.sh

echo -e "\tBuilding PON counts"
./somatic_pon_count.sh
