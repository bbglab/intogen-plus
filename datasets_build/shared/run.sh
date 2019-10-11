#!/usr/bin/env bash

set -e

if [ -z "${INTOGEN_RELEASE}" ]
then
      echo "ERROR: Define the INTOGEN_RELEASE variable"
      exit -1
fi

./ensembl_canonical_transcripts.sh
./cds.sh
./counts.sh
./somatic_pon_count.sh