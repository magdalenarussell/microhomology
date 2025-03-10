#!/bin/bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate microhomology_jax
set -eu

MOD_OUTPUT_PATH=$1
NCPU=$2
INT_MH_PARAM=$3

param_path=$(python $PWD/mh_simulation_scripts/intermediate-mh_signal_simulator_scripts/get_igor_params.py $MOD_OUTPUT_PATH)
echo "finished getting baseline igor parameters"

Rscript $PWD/mh_simulation_scripts/intermediate-mh_signal_simulator_scripts/process_data_for_intermediate-mh_signal_simulation.R $NCPU $INT_MH_PARAM
