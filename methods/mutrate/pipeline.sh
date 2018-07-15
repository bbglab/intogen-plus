#!/usr/bin/env bash

# How to use this script
# ./pipeline.sh <annotmuts_path> <genemuts_path> <mcmc_iterations> <output_folder> <cores>
# Example: ./pipeline.sh ./annotmuts.tsv ./genemuts.tsv 1000 ./output/BRCA_TCGA 30 


ANNOTMUTS=$1
GENEMUTS=$2
DNDSOUT=$3
ITERATIONS=$4
OUTPUT_FOLDER=$5
CORES=$6

if [ -d "$OUTPUT_FOLDER" ]; then
   echo "Output folder exists"
   exit 0
fi


SCRIPTS_FOLDER="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export PYTHONPATH=${SCRIPTS_FOLDER}:${PYTHONPATH}

# STEP1. Prepare files for sigfit to do signature fitting
echo "## STEP1: Preparing files for sigfit..."
mkdir -p ${OUTPUT_FOLDER}
source deactivate
source activate intogen2017_mutrate_python

if python -W ignore ${SCRIPTS_FOLDER}/mutrate.py \
                    sigfit -a ${ANNOTMUTS} \
                           -o ${OUTPUT_FOLDER}/mutational_catalogue.tsv; then
    echo "STEP1 launched successfully"
else
    echo "ERROR at STEP1"
    exit -1
fi

# STEP2. Run signature fitting with the full catalogue
echo "## STEP2: Running sigfit..."
mkdir -p ${OUTPUT_FOLDER}/sigfit_results
source deactivate
source activate intogen2017_mutrate_r
if Rscript --vanilla ${SCRIPTS_FOLDER}/cosmic_exome_fit.r -m ${OUTPUT_FOLDER}/mutational_catalogue.tsv \
                                                          -c ${SCRIPTS_FOLDER}/datasets/cosmic_exome.tsv \
                                                          -i ${ITERATIONS} \
                                                          -o ${OUTPUT_FOLDER}/sigfit_results; then
    echo "STEP2 launched successfully"
else
    echo "ERROR at STEP2"
    exit -1
fi

# STEP3. Compute mut expectation per bp as {gene : {sample : {context : expectation}}}
echo "## STEP3: Compute mut expectation per bp as {gene : {sample : {context : expectation}}}"
mkdir -p ${OUTPUT_FOLDER}/genes
source deactivate
source activate intogen2017_mutrate_python
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1
if python -W ignore ${SCRIPTS_FOLDER}/mutrate.py \
                    compute_mutrate -a ${ANNOTMUTS} \
                                    -g ${GENEMUTS} \
                                    -m ${OUTPUT_FOLDER}/mutational_catalogue.tsv \
                                    -e ${OUTPUT_FOLDER}/sigfit_results/mean_exposure.tsv \
                                    -c ${CORES} \
                                    -o ${OUTPUT_FOLDER}/genes; then
    echo "STEP3 launched successfully"
else
    # Exit on error
    echo "ERROR at STEP3"
    exit -1
fi

# STEP4. Create MAF with mutation rate
echo "## STEP4: Adding expected mutation rate to MAF table"
source deactivate
source activate intogen2017_mutrate_python
if python -W ignore ${SCRIPTS_FOLDER}/annotate_expectation.py \
                   --genes_path ${OUTPUT_FOLDER}/genes/ \
                   --maf_path ${ANNOTMUTS} \
                   --dnds_path ${DNDSOUT} \
                   --output ${OUTPUT_FOLDER}/annotmuts_expect.tsv; then
    echo "STEP4 launched successfully"
else
    # Exit on error
    echo "ERROR at STEP4"
    exit -1
fi
