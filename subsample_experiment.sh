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
L2_REG=$8
PROP=$9

Rscript $PWD/analysis_scripts/process_data_for_subsampling_experiment.R $ANNOTATION_TYPE $PARAM_GROUP $NCPU $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $MODEL_TYPE $PROP
echo "finished processing data for experiment"

python -i $PWD/jax_scripts/subsampling_experiment.py $ANNOTATION_TYPE $PARAM_GROUP $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $MODEL_TYPE $L2 $L2_REG $PROP $NCPU
echo "finished training models and making predictions"
