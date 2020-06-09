#!/bin/bash
set -e

TRANSCRIPT=$1
PFDOMAIN=$2
START=$3
END=$4
CONTAINER=$5
DATA=$6
TRANSCRIPTS_FILE=$7

range1=`singularity run -c -B ${DATA}:/data ${CONTAINER} panno --ensembl -i $TRANSCRIPT:$START | tail -n+2 | cut -f5 | cut -d'/' -f1 | sed 's/\:g\./_/g' | tr '_' '\t'`
range2=`singularity run -c -B ${DATA}:/data ${CONTAINER} panno --ensembl -i $TRANSCRIPT:$END | tail -n+2 | cut -f5 | cut -d'/' -f1 | sed 's/\:g\./_/g' | tr '_' '\t'`

gene_transcript_hugo=`cat ${TRANSCRIPTS_FILE} | grep $TRANSCRIPT`
gene=`echo -e -n "$gene_transcript_hugo" | cut -f1`
hugo=`echo -e -n "$gene_transcript_hugo" | cut -f3`
chromosome=`echo -e -n "$range1" | cut -f1 | sed 's/chr//g'`
positions=`echo -e -n "$range1\t$range2" | cut -f2,3,5,6`
min_pos=`echo -e -n "$positions" | tr '\t' '\n' | sort -n | head -n1`
max_pos=`echo -e -n "$positions" | tr '\t' '\n' | sort -n -r | head -n1`

echo -e "$chromosome\t$min_pos\t$max_pos\t-\t$TRANSCRIPT:$PFDOMAIN:$START:$END\t$TRANSCRIPT:$PFDOMAIN:$START:$END\t$hugo"
