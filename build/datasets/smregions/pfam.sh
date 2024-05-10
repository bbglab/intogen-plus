#!/usr/bin/env bash
set -e

OUTPUT=$1
BIOMART_PFAM=$2
PANNO_SCRIPT=$3
TRANSVAR=$4
TRANSVAR_DATA=$5
BIOMART_CDS=$6
CORES=$7

echo -e "CHROMOSOME\tSTART\tEND\tSTRAND\tELEMENT_ID\tSEGMENT\tSYMBOL" \
	> ${OUTPUT}

tmpdir=$(mktemp -d)

cat $BIOMART_CDS | \
	awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $11, $2}' | \
	sort -u \
	> ${tmpdir}/MANE_TRANSCRIPTS.tsv.tmp

TRANSCRIPTS=${tmpdir}/MANE_TRANSCRIPTS.tsv.tmp

if [[ "${CORES}" -eq 1 ]]
then
	zcat ${BIOMART_PFAM} > ${tmpdir}/split.0
else
	divisor=$((${CORES}-1))
	zcat ${BIOMART_PFAM} > ${tmpdir}/input
	lines=`cat ${tmpdir}/input | wc -l`
	split -l$((${lines}/${divisor})) ${tmpdir}/input \
		${tmpdir}/split. -d
fi


for f in ${tmpdir}/split*
do
	name=${f##*.}
	awk -v script="${PANNO_SCRIPT}" -v transvar="${TRANSVAR}" \
		-v transvardata="${TRANSVAR_DATA}" -v transcripts="${TRANSCRIPTS}" \
		'{system(script" "$2" "$5" "$3" "$4" "transvar" "transvardata" "transcripts)}' $f |\
		grep -v "^\s" > ${tmpdir}/${name}.tsv &
done

wait

cat ${tmpdir}/*.tsv >> ${OUTPUT}


