#!/bin/bash
set -x

# Script arguments
INPUT_FILE=$1
OUTPUT_FOLDER=$2
SIGNATURE_FILE=$3
DATA_FOLDER=$4
CORES=$5

# Enviroment variables
HYPERMUT=1000
SCRIPTS_FOLDER="/hotmaps/scripts"

DATASETS_FOLDER=$INTOGEN_DATASETS/hotmaps

# Functions

exit_on_error() {
    exit_code=$1
    if [ "${exit_code}" != "0" ]; then
        exit "${exit_code}"
    fi
}

# Preprocess
input_name=$(basename "${INPUT_FILE}")
name=${input_name%.in.maf}

TEMP_FOLDER="$OUTPUT_FOLDER/${name}.tmp"
mkdir -p "${TEMP_FOLDER}"

input_filename=${name}.maf
cp "${INPUT_FILE}" "${TEMP_FOLDER}/${input_filename}"


# TODO remove the file checks

## STEP1. Map to Structure (output: non_filtered_mupit.INPUT_FILENAME)
if [ ! -f "${TEMP_FOLDER}/non_filtered_mupit.${input_filename}" ]
then
    python $SCRIPTS_FOLDER/map_maf_to_structure.py \
        --data-dir "${TEMP_FOLDER}" \
        --match-regex "${input_filename}" \
        --output-dir "${TEMP_FOLDER}" \
        --database "${DATA_FOLDER}/mupit_database.db"
    exit_on_error $?
fi

## STEP2. Convert MAF to MUPIT (output: coverage_info.txt, INPUT_FILENAME.mupit )
if [ ! -f "${TEMP_FOLDER}/${input_filename}.mupit" ]
then
    python $SCRIPTS_FOLDER/convert_maf_to_mupit.py \
        --maf "${TEMP_FOLDER}/${input_filename}" \
        --tumor-type "${name}" \
        --no-stratify \
        -mt "${HYPERMUT}" \
        -i "${TEMP_FOLDER}" \
        --output "${TEMP_FOLDER}/${input_filename}.mupit" \
        --database "${DATA_FOLDER}/mupit_database.db"
    exit_on_error $?
fi

## STEP3. Filter hypermutated (output: mupit.INPUT_FILENAME)
if [ ! -f "${TEMP_FOLDER}/mupit.${input_filename}" ]
then
    python $SCRIPTS_FOLDER/filter_hypermutated.py \
        --raw-dir "${TEMP_FOLDER}" \
        --match-regex "${input_filename}" \
        --mut-threshold "${HYPERMUT}" \
        --sample-col 'Tumor_Sample_Barcode' \
        --data-dir "${TEMP_FOLDER}"
    exit_on_error $?
fi

## STEP4. Count mutations (input: mupit.* output: collected.INPUT_FILENAME)
if [ ! -f "${TEMP_FOLDER}/collected.${input_filename}" ]
then
    python "${SCRIPTS_FOLDER}/count_mutations.py" \
        --data-dir "${TEMP_FOLDER}"
    exit_on_error $?
fi

## STEP5. Format mutations table (input: collected.* output: mutation_tcga.INPUT_FILENAME.txt)
if [ ! -f "${TEMP_FOLDER}/mutation_tcga.${input_filename}" ]
then
    python "${SCRIPTS_FOLDER}/format_mutations_table.py" \
    	--data-dir "${TEMP_FOLDER}"
    exit_on_error $?
fi

## STEP7. Run HotMAPS (input: mutation_tcga.NAME.txt output:hotspot_INPUT_FILENAME)
if [ ! -f "${TEMP_FOLDER}/hotspot_${input_filename}" ]
then
    python "${SCRIPTS_FOLDER}/hotspot.py" \
        --log-level=INFO \
        -m "${TEMP_FOLDER}/mutation_tcga.${name}.txt" \
        -a "${DATA_FOLDER}/fully_described_pdb_info.txt" \
        -t EVERY -n 10000 -r 10.0 -c "${CORES}" \
        -o "${TEMP_FOLDER}/hotspot_${input_filename}" \
        -e "${TEMP_FOLDER}/${input_filename}.err --log=stdout" \
        -gc "${DATASETS_FOLDER}/coordinates.txt.gz" \
        -S "${SIGNATURE_FILE}" \
        --maf "${TEMP_FOLDER}/${input_filename}" \
        --database "${DATA_FOLDER}/mupit_database.db" \
        --pdb "${DATA_FOLDER}/pdb"
    exit_on_error $?
fi

## STEP8. Multiple test
if [ ! -f "${TEMP_FOLDER}/mtco_${input_filename}" ]
then
    python "${SCRIPTS_FOLDER}/multiple_testing_correction.py" \
        -i "${TEMP_FOLDER}/hotspot_${input_filename}" \
        -f min -q 0.05 \
        -m "${TEMP_FOLDER}/${input_filename}.mupit" \
        -o "${TEMP_FOLDER}/mtco_${input_filename}" \
        -s "${TEMP_FOLDER}/mtcs_${input_filename}"
    exit_on_error $?
fi

## STEP9. Find Hotspots regions gene
if [ ! -f "${TEMP_FOLDER}/hotspot_gene_${input_filename}" ]
then
    python "${SCRIPTS_FOLDER}/find_hotspot_regions_gene.py" \
        -m "${TEMP_FOLDER}/mtco_${input_filename}" \
        -a "${TEMP_FOLDER}/${input_filename}.mupit" \
        -p "${DATA_FOLDER}/fully_described_pdb_info.txt" \
        -r 10.0 -q 0.05 \
        -o "${TEMP_FOLDER}/hotspot_gene_${input_filename}" \
        --pdb "${DATA_FOLDER}/pdb"
    exit_on_error $?
fi

## STEP10. Output parser
if [ ! -f "$OUTPUT_FOLDER/${name}.out.gz" ]
then
    python "${SCRIPTS_FOLDER}/postprocess.py" \
    	"${TEMP_FOLDER}/hotspot_gene_${input_filename}" \
    	"${TEMP_FOLDER}/mtco_${input_filename}" \
    	"${OUTPUT_FOLDER}/${name}.out.gz" \
    	"${OUTPUT_FOLDER}/${name}.clusters.gz"
    exit_on_error $?
fi
