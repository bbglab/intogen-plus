#!/usr/bin/env bash
set -e

FILE=$1

name=$(basename ${FILE})
ref=${name:0:4}

wget http://hgdownload.cse.ucsc.edu/goldenPath/${ref}/liftOver/${name} \
	-O ${FILE}
# Ensure that the timestamp is updated for make
touch ${FILE}