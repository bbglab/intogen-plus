#!/usr/bin/env bash
set -e

INFO=$1
DST=$2

tmpdir=`mktemp -d`

rsync -rptz --copy-links --port=33444 ftp.rcsb.org::ftp_data/biounit/coordinates/all/*.pdb1.gz ${tmpdir}/

for it in `zcat ${INFO} | cut -f 1 | sort |uniq`
do
	file=${tmpdir}/${it}.pdb1.gz
	if [[ -f ${file} ]]
	then
		mv ${file} ${DST}
	fi
done

#rm -r ${tmpdir}
