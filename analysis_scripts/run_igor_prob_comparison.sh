#!/bin/bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate microhomology_jax
set -eu

TRAINING_DATA=$1
TEMP_DIR=$2
OUTPUT_DIR=$3
NCPU=$4

Rscript $PWD/analysis_scripts/prepare_data_for_igor.R $TRAINING_DATA 500000 $TEMP_DIR
echo "finished preparing data for igor training"

CURRENT_DIR=$PWD
cd $TEMP_DIR

python -i $PWD/igor_annotation_scripts/train/train_new_model.py $TEMP_DIR $OUTPUT_DIR 500000 alpha
echo "finished getting baseline igor parameters"

cd $CURRENT_DIR

DATA_PATH=$(Rscript $PWD/analysis_scripts/prepare_all_sequence_data.R nonproductive_v-j_trim_ligation-mh $NCPU)
echo "finished processing datasets"

python $PWD/jax_scripts/predict_twostep_EM.py $DATA_PATH igor_alpha nonproductive_v-j_trim_ligation-mh 1 2 twostep_motif_two-side-base-count-beyond_average-mh_ligation-mh True $NCPU igor_prob_experiment
echo "finished making predictions"

Rscript $PWD/plotting_scripts/manuscript_figs/analysis_igor_comparison/compare_model_probs_with_custom_igor.R nonproductive_v-j_trim_ligation-mh $NCPU 
