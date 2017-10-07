#!/bin/bash
set -x

# Script arguments
INPUT_FILE=$1
OUTPUT_FOLDER=$2
CORES=$3

# Content hotmaps.sh
# export HYPERMUT=1000
# export DATASETS_FOLDER=~/workspace/intogen/intogen-plus/datasets/hotmaps
# export DATASETS_PDB_FILE="tp53_described_pdb_info.txt"
# export MYSQL_HOST=127.0.0.1
# export MYSQL_PORT=3306
# export MYSQL_DB=mupit_modbase
# export MYSQL_USER=root
# export MYSQL_PASSWD=S82XWgESQjKJpLx3

# Preprocess
INPUT_FOLDER=$(dirname "${INPUT_FILE}")
INPUT_FILENAME=$(basename "${INPUT_FILE}")
INPUT_TUMOR_TYPE=${INPUT_FILENAME%.in.gz}
SCRIPTS_FOLDER="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/scripts
TEMP_FOLDER="$OUTPUT_FOLDER/tmp"
mkdir -p $OUTPUT_FOLDER
mkdir -p $TEMP_FOLDER

MAF_FILE=${TEMP_FOLDER}/${INPUT_TUMOR_TYPE}.maf
INPUT_FILENAME=${INPUT_TUMOR_TYPE}.maf

zcat ${INPUT_FILE} > ${MAF_FILE}

## STEP1. Map to Structure (output: non_filtered_mupit.INPUT_FILENAME)
if [ ! -f "$TEMP_FOLDER/non_filtered_mupit.$INPUT_FILENAME" ]
then
    python $SCRIPTS_FOLDER/map_maf_to_structure.py \
        --data-dir $TEMP_FOLDER \
        --match-regex $INPUT_FILENAME \
        --host $MYSQL_HOST --db $MYSQL_DB --mysql-user $MYSQL_USER --mysql-passwd S82XWgESQjKJpLx3 \
        --output-dir $TEMP_FOLDER
    rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
fi

## STEP2. Convert MAF to MUPIT (output: coverage_info.txt, INPUT_FILENAME.mupit )
if [ ! -f "$TEMP_FOLDER/$INPUT_FILENAME.mupit" ]
then
    python $SCRIPTS_FOLDER/convert_maf_to_mupit.py \
        --maf $MAF_FILE \
        -mh $MYSQL_HOST -mdb $MYSQL_DB --mysql-user $MYSQL_USER --mysql-passwd $MYSQL_PASSWD \
        --tumor-type $INPUT_TUMOR_TYPE \
        --no-stratify \
        -mt $HYPERMUT \
        -i $TEMP_FOLDER \
        --output $TEMP_FOLDER/$INPUT_FILENAME.mupit
    rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
fi

## STEP3. Filter hypermutated (output: mupit.INPUT_FILENAME)
if [ ! -f "$TEMP_FOLDER/mupit.$INPUT_FILENAME" ]
then
    python $SCRIPTS_FOLDER/filter_hypermutated.py \
        --raw-dir $TEMP_FOLDER \
        --match-regex $INPUT_FILENAME \
        --mut-threshold $HYPERMUT \
        --sample-col 'Tumor_Sample_Barcode' \
        --data-dir $TEMP_FOLDER
    rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
fi

## STEP4. Count mutations (input: mupit.* output: collected.INPUT_FILENAME)
if [ ! -f "$TEMP_FOLDER/collected.$INPUT_FILENAME" ]
then
    python $SCRIPTS_FOLDER/count_mutations.py \
        --data-dir $TEMP_FOLDER
    rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
fi

## STEP5. Format mutations table (input: collected.* output: mutation_tcga.INPUT_FILENAME.txt)
if [ ! -f "$TEMP_FOLDER/mutation_tcga.$INPUT_FILENAME" ]
then
    python $SCRIPTS_FOLDER/format_mutations_table.py --data-dir $TEMP_FOLDER
    rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
fi

## STEP7. Run HotMAPS (input: mutation_tcga.INPUT_TUMOR_TYPE.txt output:hotspot_INPUT_FILENAME)
if [ ! -f "$TEMP_FOLDER/hotspot_$INPUT_FILENAME" ]
then
    python $SCRIPTS_FOLDER/hotspot.py \
        --log-level=INFO \
        -m $TEMP_FOLDER/mutation_tcga.$INPUT_TUMOR_TYPE.txt \
        -a $DATASETS_FOLDER/$DATASETS_PDB_FILE \
        -t EVERY -n 10000 -r 10.0 -c $CORES \
        -o $TEMP_FOLDER/hotspot_$INPUT_FILENAME \
        -e $TEMP_FOLDER/$INPUT_FILENAME.err --log=stdout \
        -gc $DATASETS_FOLDER/coordinates.txt.gz \
        --maf $MAF_FILE
    rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
fi

## STEP8. Multiple test
if [ ! -f "$TEMP_FOLDER/mtco_$INPUT_FILENAME" ]
then
    python $SCRIPTS_FOLDER/multiple_testing_correction.py \
        -i $TEMP_FOLDER/hotspot_$INPUT_FILENAME \
        -f min -q 0.05 \
        -m $TEMP_FOLDER/$INPUT_FILENAME.mupit \
        -o $TEMP_FOLDER/mtco_$INPUT_FILENAME \
        -s $TEMP_FOLDER/mtcs_$INPUT_FILENAME
    rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
fi

## STEP9. Find Hotspots regions gene
if [ ! -f "$TEMP_FOLDER/hotspot_gene_$INPUT_FILENAME" ]
then
    python $SCRIPTS_FOLDER/find_hotspot_regions_gene.py \
        -m $TEMP_FOLDER/mtco_$INPUT_FILENAME \
        -a $TEMP_FOLDER/$INPUT_FILENAME.mupit \
        -p $DATASETS_FOLDER/$DATASETS_PDB_FILE \
        -r 10.0 -q 0.05 \
        -o $TEMP_FOLDER/hotspot_gene_$INPUT_FILENAME
    rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
fi

## STEP10. Output parser
if [ ! -f "$OUTPUT_FOLDER/$INPUT_TUMOR_TYPE.out.gz" ]
then
    python $SCRIPTS_FOLDER/postprocess.py $TEMP_FOLDER/hotspot_gene_$INPUT_FILENAME $TEMP_FOLDER/mtco_$INPUT_FILENAME $OUTPUT_FOLDER/$INPUT_TUMOR_TYPE.out.gz
    rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
fi