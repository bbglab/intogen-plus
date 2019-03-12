#!/bin/bash

set -xe

if [ -z "${INTOGEN_DATASETS}" ]
then
      echo "ERROR: Define the INTOGEN_DATASETS variable"
      exit -1
fi


# define the paths
path_base=$INTOGEN_DATASETS/role/
path_output=$path_base/drivers_role.tsv
# Datasets
path_base_intogen=/workspace/projects/intogen_2017/
path_cgi=$path_base/gene_MoA.tsv
path_intogen_drivers=$INTOGEN_DATASETS/drivers/vetted_drivers05.tsv
dnds_input_path=/workspace/projects/ubiquitins/run_tcga_pan/intogen_runs/20181105/dndscv/PCATLAS_WXS_PAN.out.gz


cp gene_MoA.tsv $path_cgi
python calculate_consensous_role.py --path_cgi_moa $path_cgi --path_drivers $path_intogen_drivers --path_run_dndscv $dnds_input_path \
                                    --path_output $path_output
