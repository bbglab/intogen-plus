#!/bin/bash

set -xe

if [ -z "${INTOGEN_DATASETS}" ]
then
      echo "ERROR: Define the INTOGEN_DATASETS variable"
      exit -1
fi


# define the paths
path_base=$INTOGEN_DATASETS/ptms/
path_output=$path_base/info_functional_sites.json
cp -r raw_data/ $path_base/raw_data
# parse the data
python parse_ptms.py --path_raw_data $path_base --path_output_dictionary $path_output
