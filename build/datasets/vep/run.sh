#!/usr/bin/env bash
set -e

CDS=$1
VEP_CONTAINER=$2
VEP_CACHE_DIR=$3
OUT=$4
CORES=$5

tmpdir=`mktemp -d`
echo $tmpdir

SCRIPT_FOLDER="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

python ${SCRIPT_FOLDER}/mutations.py -r ${CDS} -o ${tmpdir}/mutations.tsv

sort -k1,1 -k2n,2 -k3n,3 ${tmpdir}/mutations.tsv | uniq > ${tmpdir}/muts.tsv

if [[ "${CORES}" -eq 1 ]]
then
	mv ${tmpdir}/muts.tsv ${tmpdir}/split.0
else
	divisor=$((${CORES}-1))
	lines=`cat ${tmpdir}/muts.tsv | wc -l`
	split -l$((${lines}/${divisor})) ${tmpdir}/muts.tsv \
		${tmpdir}/split. -d
fi

for f in ${tmpdir}/split*
do
	name=${f##*.}
	singularity exec ${VEP_CONTAINER} vep -i ${f} \
		-o ${tmpdir}/split_${name}.vep.out --assembly GRCh38 \
		--no_stats --cache --offline --symbol --fork 4 \
		--protein --tab --canonical --mane --numbers \
		--no_headers \
		--dir ${VEP_CACHE_DIR} &
done

wait

# Get only MANE transcripts and most severe consequence types
cat ${tmpdir}/split_*vep.out |\
	grep "NM_" |\
	awk -F'[\t_/]' '{print($1"\t"$2"\t"$3"\t"$4"\t"$0)}' |\
	cut -f-4,8- |\
	bgzip > ${OUT}
