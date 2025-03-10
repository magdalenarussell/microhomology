#!/bin/bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate microhomology_jax
set -eu

ANNOTATION_TYPE=$1
PARAM_GROUP=$2
NCPU=$3
LEFT_MOTIF_COUNT=$4
RIGHT_MOTIF_COUNT=$5
MODEL_TYPE=$6
L2=$7
ANNOTATION_TYPE_VALIDATION=$8
PARAM_GROUP_VALIDATION=$9

python $PWD/jax_scripts/validate_twostep_model_EM.py $ANNOTATION_TYPE $PARAM_GROUP $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $MODEL_TYPE $L2 $NCPU $ANNOTATION_TYPE_VALIDATION $PARAM_GROUP_VALIDATION False
echo "finished model validation"
