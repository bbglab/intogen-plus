#!/bin/bash

set -xe

if [ -z "${INTOGEN_DATASETS}" ]
then
      echo "ERROR: Define the INTOGEN_DATASETS variable"
      exit -1
fi


# define the paths
path_base=$INTOGEN_DATASETS/stats_cohorts/
path_output=$path_base/stats_cohorts.tsv
# Datasets
path_intogen=/workspace/datasets/intogen_datasets/genomes/
path_hartwig=/workspace/datasets/hartwig/20190502/DR-024-update2/data/
path_stjude=/workspace/datasets/stjude/20180716/preprocess/
# Mapping dictionary
dict_mapping_cohorts=$path_base/dictionary_datasets.json
dict_mapping_web_names=$path_base/dictionary_long_names.json


# download the data

cp dictionary_datasets.json $dict_mapping_cohorts
cp dictionary_long_names.json $dict_mapping_web_names
python create_stats_cohorts.py --path_hartwig $path_hartwig  --path_intogen $path_intogen \
                               --path_stjude $path_stjude --dict_ttypes $dict_mapping_cohorts --path_output $path_output

