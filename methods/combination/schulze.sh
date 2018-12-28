#!/usr/bin/env bash
source activate intogen2017

# usage example
# ./new_schulze.sh ~/projects/intogen-plus/methods/schulze/output/ PCATLAS_WXS_KICH

set -x

# export SCHULZE_DATA='/workspace/projects/intogen_2017/data/latest/schulze/'

# example script arguments
# OUTPUT=~/projects/intogen-plus/methods/optimization/test_results/
# RUN_FILENAME=PCATLAS_WXS_KICH

OUTPUT=$1
RUN_FILENAME=$2

mkdir -p $OUTPUT/tmp

SCRIPTS_FOLDER="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export PYTHONPATH=$SCRIPTS_FOLDER:$PYTHONPATH

# STEP1. Read the output of intogen
echo "## STEP1"
python $SCRIPTS_FOLDER/Parser.py --input $OUTPUT/../ --cancer $RUN_FILENAME --output $OUTPUT/tmp/$RUN_FILENAME.step1

# Exit on error
rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi

# STEP2. Optimize a couple of cancer types
echo "## STEP2"
python $SCRIPTS_FOLDER/grid_optimizer.py --t_combination RANKING \
 --foutput $OUTPUT/tmp/$RUN_FILENAME.step2 \
 --input_rankings $OUTPUT/tmp/$RUN_FILENAME.step1 \
 --moutput $OUTPUT/../ \
 --cancer $RUN_FILENAME

# Exit on error
rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi

# STEP3. Generate the final combination based on the optimization
echo "## STEP3"
python $SCRIPTS_FOLDER/schulze.py \
 --input_data $OUTPUT/tmp/$RUN_FILENAME.step1 \
 --optimize_weights $OUTPUT/tmp/$RUN_FILENAME.step2 \
 --report_output $OUTPUT/tmp/$RUN_FILENAME.step3 \
 --dict_output $OUTPUT/tmp/$RUN_FILENAME.step3b

# Exit on error
rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi

# STEP4. Include the combined p-values
echo "## STEP4"
python $SCRIPTS_FOLDER/stouffer_script.py \
    --input_path $OUTPUT/tmp/$RUN_FILENAME.step1b \
    --path_rankings $OUTPUT/tmp/$RUN_FILENAME.step3 \
    --path_weights $OUTPUT/tmp/$RUN_FILENAME.step2 \
    --path_fml $OUTPUT/../oncodrivefml/$RUN_FILENAME.out.gz \
    --path_dndscv $OUTPUT/../dndscv/$RUN_FILENAME.out.gz \
    --output_path $OUTPUT/$RUN_FILENAME.stouffer.out.gz

# Exit on error
rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi

# STEP5. Generate the two tiers list at 0.01
echo "## STEP5"
python $SCRIPTS_FOLDER/create_tiers_drivers.py --threshold 0.01 \
    --column_filter QVALUE_stouffer_w \
    --input $OUTPUT/$RUN_FILENAME.stouffer.out.gz \
    --output_file $OUTPUT/$RUN_FILENAME.01.out.gz

# STEP6. Generate the two tiers list at 0.05
echo "## STEP6"
python $SCRIPTS_FOLDER/create_tiers_drivers.py --threshold 0.05 \
    --column_filter QVALUE_stouffer_w \
    --input $OUTPUT/$RUN_FILENAME.stouffer.out.gz \
    --output_file $OUTPUT/$RUN_FILENAME.05.out.gz
