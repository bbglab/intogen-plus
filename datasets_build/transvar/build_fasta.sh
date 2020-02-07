#!/bin/bash

GENOME=$1
OUTPUT=$2

BGDATA_OUT=`bgdata get datasets/genomereference/hg38`
GENOME_PATH=`echo ${BGDATA_OUT} | awk '{print $NF}'`

for file in $(ls ${GENOME_PATH}/chr*.txt); do
    basename ${file} | sed 's/chr/>/g' | sed 's/.txt//g' >> ${OUTPUT}
    fold -w61 ${file} | awk '{printf "%61s\n", $0}' >> ${OUTPUT}
done

