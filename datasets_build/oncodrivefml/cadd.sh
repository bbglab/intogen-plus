#!/bin/bash
set -e

REGIONS=$1
CADD_FILE=$2
CORES=$3
OUTPUT=$4

tmpfile=$(mktemp)

zcat $(REGIONS) | tail -n +2 |\
awk -v cadd="${CADD_URL}" '{system("tabix "cadd" "$$1":"$$2"-"$$3)}' |\
awk 'BEGIN {FS="\t";OFS = FS};{print $$1,$$2,$$3,$$4,$$6}' |\
sort --parallel=${CORES} -S 4G -k1,1 -k2,2n |\
uniq | bgzip > ${OUTPUT}
