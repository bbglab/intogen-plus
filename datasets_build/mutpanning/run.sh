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


# Create the mutpanning datasets folder
MUTPANNING_FOLDER=../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/mutpanning
mkdir -p ${MUTPANNING_FOLDER}

# Download the database
wget https://datasets.genepattern.org/data/module_support_files/MutPanning/Hg19.zip
unzip Hg19.zip
mv Hg19/ ${MUTPANNING_FOLDER}/
