#!/bin/bash

set -xe

if [ -z "${INTOGEN_DATASETS}" ]
then
      echo "ERROR: Define the INTOGEN_DATASETS variable"
      exit -1
fi


# define the paths
path_base=$INTOGEN_DATASETS/combination/
file_output=$path_base/negative_gene_set.tsv
file_output_non_expressed=$path_base/non_expressed_genes_tcga.tsv
file_olfactory_receptors=olfactory_receptors.tsv

path_expression_tcga=/workspace/datasets/TCGA_PanCanAtlas/open_version/expression/data_parsed_dataframe/


# download the data

wget https://genome.weizmann.ac.il/horde/download/genes.csv -O olfactory_receptors.tsv


# create the dataframe

python create_negative_set.py --olfactory_receptors $file_olfactory_receptors --path_expression_data $path_expression_tcga --output_total $file_output --output_non_expressed $file_output_non_expressed