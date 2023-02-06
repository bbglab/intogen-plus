#!/usr/bin/env bash
set -e

INFO=$1
DST=$2

tmpdir=`mktemp -d`

wget https://salilab.org/modbase-download/projects/genomes/H_sapiens/2013/ModBase_H_sapiens_2013_GRCh37.70.pep.all.tar.xz -P ${tmpdir}
tar -C ${tmpdir} -xf ${tmpdir}/ModBase_H_sapiens_2013_GRCh37.70.pep.all.tar.xz ModBase_H_sapiens_2013_GRCh37.70.pep.all/models/model

src=${tmpdir}/ModBase_H_sapiens_2013_GRCh37.70.pep.all/models/model

for prot in `zcat ${INFO} | cut -f 1 | grep "^ENSP" | sort |uniq`
do
	for file in ${src}/${prot}*
	do
		filename="$(basename ${file})"
		if [[ ${file} =~ \.xz ]]; then
			xz -cd ${file} | gzip > ${DST}/${filename/\.xz/\.gz}
		elif [[ ${file} =~ \.gz ]]; then
			mv ${file} ${DST}/
		else
			gzip -c ${file} > ${DST}/${filename}.gz
		fi
	done
done

#rm -r ${tmpdir}