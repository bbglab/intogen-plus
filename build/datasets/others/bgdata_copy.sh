#!/usr/bin/env bash
set -e
FOLDER=$1

bgdata_out=`bgdata get intogen/pon/hmf`
data_path=`echo ${bgdata_out} | awk '{print $NF}'`
cp ${data_path} ${FOLDER}/