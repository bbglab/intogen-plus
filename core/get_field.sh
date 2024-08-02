#!/bin/bash
set -e

input=$1
field=$2

C=1
header=$( zcat "${input}" | head -n 1 )
for i in ${header}
do
	if [[ "$i" == "${field}" ]]
	then
		break
	else
		C=$(( C + 1 ))
	fi
done

zcat "${input}" | sed -n '2p' | cut -f$C | xargs printf