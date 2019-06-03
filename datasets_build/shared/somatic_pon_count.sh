#!/usr/bin/env bash

set -e

if [ -z "${INTOGEN_RELEASE}" ]
then
      echo "ERROR: Define the INTOGEN_RELEASE variable"
      exit -1
fi


INPUT_URL="https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd/download?path=%2FHMF-Pipeline-Resources&files=SOMATIC_PON.vcf.gz"
OUTPUT_FILE="../../datasets/${INTOGEN_RELEASE}/shared/somatic_pon_count_filtered.tsv.gz"

python somatic_pon_counts.py -u ${INPUT_URL} -o ${OUTPUT_FILE}
