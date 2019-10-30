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


# define the paths
INTOGEN_DATASETS=../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}
path_base=${INTOGEN_DATASETS}/combination/
file_output=${path_base}/negative_gene_set.tsv
file_output_non_expressed=${path_base}/non_expressed_genes_tcga.tsv
file_olfactory_receptors=${path_base}/olfactory_receptors.tsv

mkdir -p ${path_base}

# path_expression_tcga=/workspace/datasets/TCGA_PanCanAtlas/open_version/expression/data_parsed_dataframe/
bgdata get --force intogen/expression/tcga_pancanatlas

# download the data
wget https://genome.weizmann.ac.il/horde/download/genes.csv -O ${file_olfactory_receptors}

# create the dataframe
python create_negative_set.py \
  --olfactory_receptors ${file_olfactory_receptors} \
  --output_total ${file_output} \
  --output_non_expressed ${file_output_non_expressed}
