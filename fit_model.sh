#!/bin/bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate mechanistic-trimming_jax
set -eu

ANNOTATION_TYPE=$1
PARAM_GROUP=$2
NCPU=$3
LEFT_MOTIF_COUNT=$4
RIGHT_MOTIF_COUNT=$5
MODEL_TYPE=$6
L2=$7

Rscript $PWD/mechanistic-trimming/scripts/process_data_for_model_fitting.R $ANNOTATION_TYPE $PARAM_GROUP $NCPU $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $MODEL_TYPE
echo "finished processing data for model fitting"

python -i $PWD/mechanistic-trimming/jax_scripts/fit_model.py $ANNOTATION_TYPE $PARAM_GROUP $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $MODEL_TYPE $L2 $NCPU
echo "finished training model and making predictions"
