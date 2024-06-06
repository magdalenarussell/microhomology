#!/bin/bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate microhomology_jax
set -eu

MOD_OUTPUT_PATH=$1
PARAM_GROUP=$2
NCPU=$3
NP_COND=$4

param_path=$(python $PWD/analysis_scripts/trim_lig_config_mh_simulator/get_igor_params.py $MOD_OUTPUT_PATH)
echo "finished getting baseline igor parameters"

Rscript $PWD/analysis_scripts/trim_lig_config_mh_simulator/process_data_for_ligation-mh_trim_signal_simulation.R $PARAM_GROUP $NCPU $NP_COND
