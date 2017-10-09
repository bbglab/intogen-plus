#!/usr/bin/env bash

source activate intogen2017
# export SCHULZE_DATA = '/home/jordeu/workspace/intogen/intogen-plus/datasets/schulze'
# export OUTPUT=/home/jordeu/workspace/intogen/intogen-plus/output
# export RUN_FILENAME=KICH.out.gz


# STEP1. Read the output of intogen
python Parser.py --input $OUTPUT --cancer KICH --output $OUTPUT/schulze/$RUN_FILENAME.step1

# STEP2. Optimize a couple of cancer types
python optimizer.py --seeds 1 --niter 1 --epsilon 0.1 --t_combination RANKING \
 --foutput $OUTPUT/schulze/$RUN_FILENAME.step2 \
 --input_rankings $OUTPUT/schulze/$RUN_FILENAME.step1 \
 --discarded_methods $SCHULZE_DATA/discarted_analyses.txt

# STEP3. Generate the final combination based on the optimization
python schulze.py --name Combination_Ranking_Optimized --type_input ranking \
 --input_data $OUTPUT/schulze/$RUN_FILENAME.step1 \
 --optimize_weights $OUTPUT/schulze/$RUN_FILENAME.step2 \
 --directory_output $OUTPUT/schulze \
 --dict_output $OUTPUT/dict_ranking_optimized.pickle

# STEP4. Include the combined p-values
python stouffer_script.py \
    --tumor_type BRCA \
    --input_path $OUTPUT/schulze/BRCA.tsv \
    --output_path $OUTPUT/schulze/BRCA_combined.tsv \
    --path_ranking $OUTPUT/schulze/BRCA_ranking.tsv \
    --path_weights $OUTPUT/schulze/BRCA_weights.tsv \
    --path_fml $OUTPUT/oncodrivefml

# STEP5. Generate the two tiers list
python create_tiers_drivers.py --threshold 0.01 --column_filter QVALUE_stouffer_w \
    --input_dir $OUTPUT/schulze \
    --output_file /workspace/projects/intogen/test/output_tiers.tsv \
    --pattern _combined.tsv

