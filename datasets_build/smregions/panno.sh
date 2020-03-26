#!/bin/bash

TRANSCRIPT=$1
PFDOMAIN=$2
START=$3
END=$4
CONTAINER=$5
DATA=$6
TRANSCRIPTS_FILE=$7


INTOGEN_DATASETS="../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}"

RANGE1=`singularity run -c -B ${DATA}:/data ${CONTAINER} panno --ensembl -i $TRANSCRIPT:$START | tail -n+2 | cut -f5 | cut -d'/' -f1 | sed 's/\:g\./_/g' | tr '_' '\t'`
RANGE2=`singularity run -c -B ${DATA}:/data ${CONTAINER} panno --ensembl -i $TRANSCRIPT:$END | tail -n+2 | cut -f5 | cut -d'/' -f1 | sed 's/\:g\./_/g' | tr '_' '\t'`

GENE_TRANSCRIPT_HUGO=`cat ${TRANSCRIPTS_FILE} | grep $TRANSCRIPT`
GENE=`echo -e -n "$GENE_TRANSCRIPT_HUGO" | cut -f1`
HUGO=`echo -e -n "$GENE_TRANSCRIPT_HUGO" | cut -f3`
CHROMOSOME=`echo -e -n "$RANGE1" | cut -f1 | sed 's/chr//g'`
POSITIONS=`echo -e -n "$RANGE1\t$RANGE2" | cut -f2,3,5,6`
MIN_POS=`echo -e -n "$POSITIONS" | tr '\t' '\n' | sort -n | head -n1`
MAX_POS=`echo -e -n "$POSITIONS" | tr '\t' '\n' | sort -n -r | head -n1`

echo -e "$CHROMOSOME\t$MIN_POS\t$MAX_POS\t-\t$TRANSCRIPT:$PFDOMAIN:$START:$END\t$TRANSCRIPT:$PFDOMAIN:$START:$END\t$HUGO"
