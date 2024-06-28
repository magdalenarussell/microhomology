#!/bin/bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate microhomology_jax
set -eu

TRAINING_DATA=$1
TEMP_DIR=$2
SAMPLE_COUNT=$3

IGOR_DATA=$(Rscript $PWD/analysis_scripts/igor_infer/prepare_data_for_igor.R $TRAINING_DATA $SAMPLE_COUNT $TEMP_DIR)
echo "finished preparing data for igor training"

CURRENT_DIR=$PWD
cd $TEMP_DIR

python $CURRENT_DIR/analysis_scripts/igor_infer/compile_data_and_run.py $IGOR_DATA $CURRENT_DIR
echo "finished inferring IGoR model"

cd $CURRENT_DIR
