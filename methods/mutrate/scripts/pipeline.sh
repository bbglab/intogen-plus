#!/usr/bin/env bash

# How to use this script
# ./pipeline.sh <annotmuts_path> <genemuts_path> <mcmc_iterations> <job_id> <cores>
# Example: ./pipeline.sh ./annotmuts.tsv ./genemuts.tsv 1000 BRCA_TCGA 30

source activate signatures_env

PATH_BASE="/home/fmuinos/projects/intogen-plus/methods/mutrate"
PATH_OUT="/home/fmuinos/projects/intogen-plus/methods/mutrate/pipeline_results"
PATH_TMP=$(mktemp -d -p ${PATH_BASE})

ANNOTMUTS=$1
GENEMUTS=$2
ITERATIONS=$3
JOB_ID=$4
CORES=$5

export PYTHONPATH=${PATH_BASE}:${PYTHONPATH}

# STEP1. Prepare files for sigfit to do signature fitting
echo "## STEP1: Preparing files for sigfit..."
mkdir -p ${PATH_BASE}/ppln_results/${JOB_ID}
if python -W ignore ${PATH_BASE}/scripts/mutrate.py \
                    sigfit -a ${ANNOTMUTS} \
                           -o ${PATH_BASE}/ppln_results/${JOB_ID}/mutational_catalogue.tsv; then
    echo "STEP1 launched successfully"
else
    # Exit on error
    echo "$ Error in JOB ${JOB_ID} at STEP1" >> ${PATH_BASE}/failed.txt;
    rm -r ${PATH_TMP}
    exit 1  # terminate and indicate error
fi

# STEP2. Run signature fitting with the full catalogue
echo "## STEP2: Running sigfit..."
mkdir -p ${PATH_BASE}/ppln_results/${JOB_ID}/sigfit_results
if Rscript --vanilla ${PATH_BASE}/scripts/cosmic_exome_fit.r -m ${PATH_BASE}/ppln_results/${JOB_ID}/mutational_catalogue.tsv \
                                                             -c ${PATH_BASE}/cosmic_exome.tsv \
                                                             -i ${ITERATIONS} \
                                                             -o ${PATH_BASE}/ppln_results/${JOB_ID}/sigfit_results; then
    echo "STEP2 launched successfully"
else
    # Exit on error
    echo "$ Error in JOB ${JOB_ID} at STEP2" >> ${PATH_BASE}/failed.txt;
    rm -r ${PATH_TMP}
    exit 1  # terminate and indicate error
fi

# STEP3. Compute mut expectation per bp as {gene : {sample : {context : expectation}}}
echo "## STEP3: Compute mut expectation per bp as {gene : {sample : {context : expectation}}}"
if python -W ignore ${PATH_BASE}/scripts/mutrate.py \
                    compute_mutrate -a ${ANNOTMUTS} \
                                    -g ${GENEMUTS} \
                                    -m ${PATH_BASE}/ppln_results/${JOB_ID}/mutational_catalogue.tsv \
                                    -e ${PATH_BASE}/ppln_results/${JOB_ID}/sigfit_results/mean_exposure.tsv \
                                    -c ${CORES} \
                                    -o ${PATH_BASE}/ppln_results/${JOB_ID}/${JOB_ID}.out.pickle.gz; then
    echo "STEP3 launched successfully"
else
    # Exit on error
    echo "$ Error in JOB ${JOB_ID} at STEP3" >> ${PATH_BASE}/failed.txt;
    rm -r ${PATH_TMP}
    exit 1 # terminate and indicate error
fi

rm -r ${PATH_TMP}
