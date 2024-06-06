#!/bin/bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate microhomology_jax
set -eu

MOD_OUTPUT_PATH=$1
PARAM_GROUP=$2
NCPU=$3
TRIM_SAMPLING_TYPE=$4
LIGATION_MH_PARAM=$5
INT_MH_PARAM=$6
NP_COND=$7

param_path=$(python $PWD/mh_simulation_scripts/ligation-mh_signal_simulator/get_igor_params.py $MOD_OUTPUT_PATH)
echo "finished getting baseline igor parameters"

Rscript $PWD/mh_simulation_scripts/ligation-mh_signal_simulator/process_data_for_ligation-mh_trim_signal_simulation.R $PARAM_GROUP $NCPU $TRIM_SAMPLING_TYPE $LIGATION_MH_PARAM $INT_MH_PARAM $NP_COND
