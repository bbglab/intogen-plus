#!/usr/bin/env bash

FOLDER=$1

tmpdir=`mktemp -d`

wget https://datasets.genepattern.org/data/module_support_files/MutPanning/Hg19.zip \
	-O ${tmpdir}/mutpanning.zip
unzip -d ${tmpdir} ${tmpdir}/mutpanning.zip
mv ${tmpdir}/Hg19 ${FOLDER}/
chmod -R g+r ${FOLDER}/Hg19
