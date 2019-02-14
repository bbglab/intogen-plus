#!/bin/bash

set -xe

if [ -z "${INTOGEN_DATASETS}" ]
then
      echo "ERROR: Define the INTOGEN_DATASETS variable"
      exit -1
fi


# Pipeline to compute the driver list from the output of all datasets of intogen

set -e

# Pipeline to calculate features of degrons
# ---------------------------

# Help prompt

if [ "$1" == "--help" ]; then
	echo "Usage: bash get_driver_list `basename $0`"
	echo ""
	echo "Pipeline to compute the driver list from the output of intogen. "
	echo "Check the paths in the current script to edit the source files. "
	exit 0
fi

# Base paths

base_intogen=/workspace/projects/intogen_2017/runs/20180903/
base_hartwig=/workspace/projects/hartwig/intogen/runs/20180529_20180903/
base_stjude=/workspace/projects/stjude/intogen/runs/20180919/

# deconstruct for vetting

intogen_decons=${base_intogen}/deconstructsig/
hartwig_decons=${base_hartwig}/deconstructsig/
stjude_decons=${base_stjude}/deconstructsig/

# drivers to compute the final list

intogen_paths=${base_intogen}/combination/
hartwig_paths=${base_hartwig}/combination/
stjude_paths=${base_stjude}/combination/

# Path with the information of the input cohorts

info_cohorts=${INTOGEN_DATASETS}/stats_cohorts/stats_cohorts.tsv

# Path with the CGC information

cgc_info=${INTOGEN_DATASETS}/combination/cgc/cancer_gene_census_parsed.tsv

# Path scripts

path_script_vetting=prepare_vetting_files.py
path_script_drivers=get_driver_list.py

# Path output

path_output=${INTOGEN_DATASETS}/drivers/
output_vetting=${path_output}/information_vetting_genes.tsv.gz
output_intogen=${path_output}/drivers_intogen.05.tsv

# Download exact

echo "Downloading ExAC"
wget https://storage.googleapis.com/gnomad-public/release/2.1/ht/constraint/constraint.txt.bgz
mv constraint.txt.bgz $INTOGEN_DATASETS/drivers/constraint.txt.gz

exit

exact_file=$INTOGEN_DATASETS/drivers/constrains.txt.gz

expression_file_tcga=$INTOGEN_DATASETS/combination/non_expressed_genes_tcga.tsv

# Prepare the olfactory receptors

or_path=$INTOGEN_DATASETS/drivers/olfactory_receptors.tsv

cp olfactory_receptors.tsv $or_path

# Run the vetting preparison

echo "Preparing vetting DataFrame. This might take a while..."
echo $exact_file

if [ ! -f ${output_vetting}  ]; then

     python ${path_script_vetting}  -i ${intogen_decons} -w ${hartwig_decons} -s ${stjude_decons}  -c ${info_cohorts} -o ${output_vetting} \
    -t $expression_file_tcga -e $exact_file -r $or_path

else
    echo "File ${output_vetting} already exists! Skipping preparison of vetting files. If you what to recompute it, remove the output file"
fi



echo "Now running the driver discovery from the output of intOGen and the vetting information...."

if [ ! -f ${output_intogen}  ]; then


    python ${path_script_drivers}  -i ${intogen_paths} -w ${hartwig_paths} -s ${stjude_paths}  -c ${info_cohorts} -o ${path_output} -t 05 -v ${output_vetting} -g ${cgc_info}

else
    echo "File ${output_intogen} already exists! Skipping driver list discovery. If you what to recompute it, rename the output file"
fi


echo ""

echo "Done..."

exit 0