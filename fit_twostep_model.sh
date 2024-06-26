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
ITER=${8:-"NA"}

if [ $ITER = "NA" ]; then
    python -i $PWD/jax_scripts/fit_twostep_model.py $ANNOTATION_TYPE $PARAM_GROUP $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $MODEL_TYPE $L2 $NCPU
    echo "finished training model and making predictions"
else
    python -i $PWD/jax_scripts/fit_twostep_model.py $ANNOTATION_TYPE $PARAM_GROUP $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $MODEL_TYPE $L2 $NCPU $ITER
    echo "finished training model and making predictions"
fi
