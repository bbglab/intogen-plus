#!/bin/bash

set -xe

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

GENOME=$(mktemp --suffix=.fa)
for f in $(ls ../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/bgdata/datasets/genomereference/hg38-20161209/chr*.txt); do
    basename $f | sed 's/chr/>/g' | sed 's/.txt//g' >> ${GENOME}
    fold -w61 $f | awk '{printf "%61s\n", $0}' >> ${GENOME}
done

echo "library(dndscv); buildref(\"../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/shared/cds_biomart.tsv\", \"${GENOME}\", outfile = \"../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/RefCDS.rda\")" | singularity exec ../../containers/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/dndscv.simg R --no-save

rm ${GENOME}