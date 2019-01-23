#!/bin/bash

set -e

# define the paths
path_base="/workspace/projects/intogen_2017/"
path_data_intogen=$path_base/data/2019_01_18/stats_cohorts.tsv
path_cgc=$path_base/data/2019_01_21/cgc
file_cgc=$path_cgc/cancer_gene_census.csv
dict_mapping_cgc=$path_base/data/latest/mapping_cgc_ttypes.json
dict_mapping_cgc_intogen=$path_base/data/latest/mapping_cgc_ttypes_intogen.json
path_output=$path_base/data/latest/

# download the data

echo "python download_cgc.py --download path_cgc"

#python download_cgc.py --download $path_cgc

echo "cgc file store in "$file_cgc

# parse the data

echo "Parsing the dataframe"
echo "python parse_cgc.py --path_cgc_original $file_cgc --path_cohorts $path_data_intogen \
                    --dict_mapping_cgc $dict_mapping_cgc --dict_mapping_cgc_intogen $dict_mapping_cgc_intogen \
                    --path_output $path_output --debug"

python parse_cgc.py --path_cgc_original $file_cgc --path_cohorts $path_data_intogen \
                    --dict_mapping_cgc $dict_mapping_cgc --dict_mapping_cgc_intogen $dict_mapping_cgc_intogen \
                    --path_output $path_output --debug
