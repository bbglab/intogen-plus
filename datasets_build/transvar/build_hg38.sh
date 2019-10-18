#!/bin/bash

set -xe

ENSEMBL_RELEASE=${INTOGEN_VEP/vep/}
ENSEMBL_REFERENCE=${INTOGEN_GENOME/hg/GRCh}
ENSEMBL_REFERENCE=${ENSEMBL_REFERENCE/GRCh19/GRCh37}

GENOME="Homo_sapiens.${ENSEMBL_REFERENCE}.fa"
echo -e "\tBuild genome fasta file"
LAST_VERSION=`cat ../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/bgdata/datasets/genomereference/hg38.master`
for f in $(ls ../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/bgdata/datasets/genomereference/hg38-${LAST_VERSION}/chr*.txt); do
    basename $f | sed 's/chr/>/g' | sed 's/.txt//g' >> ${GENOME}
    fold -w61 $f | awk '{printf "%61s\n", $0}' >> ${GENOME}
done

# Index fasta file
TRANSVAR_DATA="../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/transvar"
mkdir -p $TRANSVAR_DATA
TRANSVAR="singularity run -B ${TRANSVAR_DATA}:/data ../../containers/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/transvar.simg"

echo -e "\tIndexing fasta file"
$TRANSVAR index --reference ${GENOME}

echo -e "\tConfigure genome reference"
mv ${GENOME}* $TRANSVAR_DATA
$TRANSVAR config -k reference -v /data/${GENOME} --refversion ${ENSEMBL_REFERENCE}

echo -e "\tDownloading ENSEMBL GTF"
ENSEMBL_GTF="Homo_sapiens.${ENSEMBL_REFERENCE}.${ENSEMBL_RELEASE}.gtf.gz"
wget "ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/homo_sapiens/${ENSEMBL_GTF}"

echo -e "\tIndexing ENSEMBL"
$TRANSVAR index --ensembl ${ENSEMBL_GTF}

echo -e "\tConfigure ENSEMBL"
mv ${ENSEMBL_GTF}.transvardb* $TRANSVAR_DATA/
rm ${ENSEMBL_GTF}*
$TRANSVAR config -k ensembl -v /data/${ENSEMBL_GTF}.transvardb --refversion ${ENSEMBL_REFERENCE}
