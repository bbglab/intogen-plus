#!/usr/bin/env bash

set -xe

if [ -z "${INTOGEN_DATASETS}" ]
then
      echo "ERROR: Define the INTOGEN_RELEASE variable"
      exit -1
fi


ANNOTMUTS=$1
GENEMUTS=$2
WEIGHTS=$3
CORES=$4
OUTPUT=$5


/usr/bin/python3 /mutrate/compute_mutrate.py compute-mutrate \
                        --annotmuts ${ANNOTMUTS} \
                        --genemuts ${GENEMUTS} \
                        --weights ${WEIGHTS} \
                        --cores ${CORES} \
                        --output ${OUTPUT}



FOLDER=`dirname ${WEIGHTS}`
NAME=`basename ${WEIGHTS} | cut -d '.' -f1`
ASSIGNMENT_PATH=${FOLDER}/${NAME}.signature_likelihood

/usr/bin/python3 /mutrate/signature_assignment.py \
                        --input_file ${WEIGHTS} \
                        --assignment_path ${ASSIGNMENT_PATH}