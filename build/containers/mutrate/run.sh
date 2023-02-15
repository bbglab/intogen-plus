#!/usr/bin/bash

set -xe

if [ -z "${INTOGEN_DATASETS}" ]
then
      echo "ERROR: Define the INTOGEN_RELEASE variable"
      exit -1
fi


ANNOTMUTS=$1
SCOPE=$2
OUTPUT=$3

python3 /mutrate/compute_mutrate.py mutrate \
                        --annotmuts ${ANNOTMUTS} \
                        --scope ${SCOPE} \
                        --outfolder ${OUTPUT}
