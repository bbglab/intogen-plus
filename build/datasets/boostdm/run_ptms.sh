#!/bin/bash

set -xe

if [ -z "${INTOGEN_DATASETS}" ]
then
      echo "ERROR: Define the INTOGEN_DATASETS variable"
      exit -1
fi


# define the paths
DEST=$1
SOURCE=$2

path_output=${DEST}/info_functional_sites.json

unzip -o ${SOURCE}/phosphosite.zip -d ${DEST}

# parse the data
python ${SOURCE}/parse_ptms.py --path_raw_data $DEST/raw_data --path_output_dictionary $path_output