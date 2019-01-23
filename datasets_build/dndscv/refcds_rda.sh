#!/bin/bash

set -xe

if [ -z "${INTOGEN_RELEASE}" ]
then
      echo "ERROR: Define the INTOGEN_RELEASE variable"
      exit -1
fi

GENOME=$(mktemp --suffix=.fa)
for f in $(ls ../../datasets/latest/bgdata/datasets/genomereference/hg38-20161209/chr*.txt); do 
    basename $f | sed 's/chr/>/g' | sed 's/.txt//g' >> ${GENOME}
    fold -w61 $f | awk '{printf "%61s\n", $0}' >> ${GENOME}
done

echo "library(dndscv); buildref(\"../../datasets/${INTOGEN_RELEASE}/shared/cds_biomart.tsv\", \"${GENOME}\", outfile = \"../../datasets/${INTOGEN_RELEASE}/RefCDS.rda\")" | singularity exec ../../containers/${INTOGEN_RELEASE}/dndscv.simg R --no-save

rm ${GENOME}