#!/bin/bash

GENOME=$1
OUTPUT=$2

bgdata_out=`bgdata get datasets/genomereference/${GENOME}`
genome_path=`echo ${bgdata_out} | awk '{print $NF}'`

for file in $(ls ${genome_path}/chr*.txt); do
    basename ${file} | sed 's/chr/>/g' | sed 's/.txt//g' >> ${OUTPUT}
    fold -w61 ${file} | awk '{printf "%61s\n", $0}' >> ${OUTPUT}
done

