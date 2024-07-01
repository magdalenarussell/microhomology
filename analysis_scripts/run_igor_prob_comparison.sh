#!/bin/bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate microhomology_jax
set -eu

TRAINING_DATA_PATH=$1
OUTPUT_PATH=$2
NCPU=$3

bash $PWD/analysis_scriptes/igor_infer/process_data_and_run_igor.sh $TRAINING_DATA_PATH $OUTPUT_PATH all

IGOR_MARG_FILE_PATH="${OUTPUT_PATH}/igor_sampled_sequences_all_igor_infer/batchname_inference/final_marginals.txt"
IGOR_PARAM_FILE_PATH="${OUTPUT_PATH}/igor_sampled_sequences_all_igor_infer/batchname_inference/final_parms.txt"

param_path=$(python $PWD/analysis_scripts/get_igor_params.py $IGOR_MARG_FILE_PATH $IGOR_PARAM_FILE_PATH $OUTPUT_PATH)
echo "finished getting baseline igor parameters"

DATA_PATH=$(Rscript $PWD/analysis_scripts/prepare_all_sequence_data.R nonproductive_v-j_trim_ligation-mh $NCPU)
echo "finished processing datasets"

python $PWD/jax_scripts/predict_twostep_EM.py $DATA_PATH igor_alpha nonproductive_v-j_trim_ligation-mh 1 2 twostep_motif_two-side-base-count-beyond_average-mh_ligation-mh True $NCPU igor_prob_experiment
echo "finished making predictions"

Rscript $PWD/plotting_scripts/manuscript_figs/analysis_igor_comparison/compare_model_probs_with_igor.R 
