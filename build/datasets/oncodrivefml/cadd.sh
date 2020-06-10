#!/bin/bash
set -e

REGIONS=$1
CADD_FILE=$2
CORES=$3
OUTPUT=$4

tmpdir=$(mktemp -d)

if [[ "${CORES}" -eq 1 ]]
then
	zcat ${REGIONS} | tail -n +2 > ${tmpdir}/split.0
else
	divisor=$((${CORES}-1))
	zcat ${REGIONS} | tail -n +2 > ${tmpdir}/input
	lines=`cat ${tmpdir}/input | wc -l`
	split -l$((${lines}/${divisor})) ${tmpdir}/input \
		${tmpdir}/split. -d
fi

for f in ${tmpdir}/split*
do
	name=${f##*.}
	awk -v cadd="${CADD_FILE}" '{system("tabix "cadd" "$1":"$2"-"$3)}' $f |\
		awk 'BEGIN {FS="\t";OFS = FS};{print $1,$2,$3,$4,$6}' |\
		> ${tmpdir}/${name}.tsv &
done

wait

cat ${tmpdir}/*.tsv |\
	sort --parallel=${CORES} -S 4G -k1,1 -k2,2n |\
	uniq | bgzip > ${OUTPUT}
