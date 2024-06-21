#!/bin/bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate microhomology_jax
set -eu

MOD_OUTPUT_PATH=$1
NCPU=$3

param_path=$(python $PWD/mh_simulation_scripts/ligation-mh_signal_simulator/get_igor_params.py $MOD_OUTPUT_PATH)
echo "finished getting baseline igor parameters"

DATA_PATH=$(Rscript $PWD/analysis_scripts/prepare_all_sequence_data.R nonproductive_v-j_trim_ligation-mh $NCPU)
echo "finished processing datasets"

python $PWD/jax_scripts/predict_twostep_EM.py $DATA_PATH igor_alpha nonproductive_v-j_trim_ligation-mh 1 2 twostep_motif_two-side-base-count-beyond_average-mh_ligation-mh True $NCPU igor_prob_experiment
echo "finished making predictions"

Rscript $PWD/plotting_scripts/manuscript_figs/analysis_igor_comparison/compare_model_probs_with_igor.R nonproductive_v-j_trim_ligation-mh $NCPU 
