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
	echo "Usage: bash boostDM_features `basename $0`"
	echo ""
	echo "Pipeline to aggregate features output from all intogen cohort "
	echo "Check the paths in the current script to edit the source files. "
	exit 0
fi

# Base paths

base_intogen=/workspace/projects/intogen_2017/runs/20190325/
base_hartwig=/workspace/projects/hartwig/intogen/runs/20181011_20190325/
base_stjude=/workspace/projects/stjude/intogen/runs/20190325/

# clustl for vetting

intogen_clustl=${base_intogen}/oncodriveclustl/*.clusters.gz
hartwig_clustl=${base_hartwig}/oncodriveclustl/*.clusters.gz
stjude_clustl=${base_stjude}/oncodriveclustl/*.clusters.gz

path_intogen_clustl=$INTOGEN_DATASETS/boostDM/oncodriveclustl/

if [ ! -d ${path_intogen_clustl}  ]; then
    mkdir $path_intogen_clustl
fi

ln -s $intogen_clustl $path_intogen_clustl
ln -s $hartwig_clustl $path_intogen_clustl
ln -s $stjude_clustl $path_intogen_clustl

# hotmaps to compute the final list

intogen_hotmaps=${base_intogen}/hotmaps/*.clusters.gz
hartwig_hotmaps=${base_hartwig}/hotmaps/*.clusters.gz
stjude_hotmaps=${base_stjude}/hotmaps/*.clusters.gz

path_intogen_hotmaps=$INTOGEN_DATASETS/boostDM/hotmaps/

if [ ! -d ${path_intogen_hotmaps}  ]; then
    mkdir $path_intogen_hotmaps
fi

ln -s $intogen_hotmaps $path_intogen_hotmaps
ln -s $hartwig_hotmaps $path_intogen_hotmaps
ln -s $stjude_hotmaps $path_intogen_hotmaps

# smregions to compute the final list

intogen_smregions=${base_intogen}/smregions/*.out.gz
hartwig_smregions=${base_hartwig}/smregions/*.out.gz
stjude_smregions=${base_stjude}/smregions/*.out.gz

path_intogen_smregions=$INTOGEN_DATASETS/boostDM/smregions/

if [ ! -d ${path_intogen_smregions}  ]; then
    mkdir $path_intogen_smregions
fi

ln -s $intogen_smregions $path_intogen_smregions
ln -s $hartwig_smregions $path_intogen_smregions
ln -s $stjude_smregions $path_intogen_smregions



echo ""

echo "Done..."

exit 0