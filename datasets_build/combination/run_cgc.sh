#!/bin/bash

set -xe

if [ -z "${INTOGEN_RELEASE}" ]
then
      echo "ERROR: Define the INTOGEN_RELEASE variable"
      exit -1
fi

if [ -z "${INTOGEN_VEP}" ]
then
      echo "ERROR: Define the INTOGEN_VEP variable"
      exit -1
fi

if [ -z "${INTOGEN_GENOME}" ]
then
      echo "ERROR: Define the INTOGEN_GENOME variable"
      exit -1
fi

INTOGEN_DATASETS=../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}

# define the paths
path_base=${INTOGEN_DATASETS}/combination/
path_cgc=${INTOGEN_DATASETS}/combination/cgc/
file_cgc=${path_base}/cgc/cancer_gene_census.csv
dict_mapping_cgc=${path_base}/cgc/mapping_cgc_ttypes.json
dict_mapping_cgc_intogen=${path_base}/cgc/mapping_cgc_ttypes_intogen.json


# download the data
echo "python download_cgc.py --download ${path_base}"
python download_cgc.py --download ${path_base}
echo "cgc file store in "${file_cgc}

# parse the data
cp mapping_cgc_ttypes.json ${INTOGEN_DATASETS}/combination/cgc
cp mapping_cgc_ttypes_intogen.json ${INTOGEN_DATASETS}/combination/cgc
# cp negative_gene_set.tsv ${INTOGEN_DATASETS}/combination

echo "Parsing the dataframe"
echo "python parse_cgc.py \\
  --path_cgc_original ${file_cgc} \\
  --dict_mapping_cgc ${dict_mapping_cgc} \\
  --dict_mapping_cgc_intogen ${dict_mapping_cgc_intogen} \\
  --path_output ${path_output}
  --debug"

python parse_cgc.py \
  --path_cgc_original ${file_cgc} \
  --dict_mapping_cgc ${dict_mapping_cgc} \
  --dict_mapping_cgc_intogen ${dict_mapping_cgc_intogen} \
  --path_output ${path_cgc} \
  --debug
