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

echo -e "\tDownloading liftover chain files"
PREPROCESS=../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/preprocess
mkdir -p ${PREPROCESS}
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz -O ${PREPROCESS}/hg19ToHg38.over.chain.gz
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg38ToHg19.over.chain.gz -O ${PREPROCESS}/hg38ToHg19.over.chain.gz

echo -e "\tCopy pileup mappability files"
cp hg19_100bp.coverage.regions.gz ${PREPROCESS}/.
cp hg38_100bp.coverage.regions.gz ${PREPROCESS}/.
