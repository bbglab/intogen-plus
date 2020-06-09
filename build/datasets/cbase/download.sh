#!/usr/bin/env bash

FOLDER=$1

tmpdir=`mktemp -d`

wget -c http://genetics.bwh.harvard.edu/cbase/CBaSE_v1.1.zip \
		-O ${tmpdir}/cbase.zip
unzip -d ${tmpdir} ${tmpdir}/cbase.zip
mv ${tmpdir}/CBaSE_v1.1/Auxiliary/*.gz ${FOLDER}/
mv ${tmpdir}/CBaSE_v1.1/Auxiliary/*.txt ${FOLDER}/