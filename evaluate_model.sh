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
CALCULATE_EXPECTED_LOSS=$8

python -i $PWD/mechanistic-trimming/jax_scripts/validate_model.py $ANNOTATION_TYPE $PARAM_GROUP $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $MODEL_TYPE $L2 $NCPU $ANNOTATION_TYPE $CALCULATE_EXPECTED_LOSS
echo "finished model validation"
