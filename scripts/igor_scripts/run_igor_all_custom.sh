#!/bin/bash

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate mechanistic-trimming_py 
echo "set environ"
set -eu

RAW_FILE=$1
TEMP_DIR=$2
OUTPUT_DIR=$3
MODEL_PARAMS=$4
MODEL_MARG=$5
ANNOTATION_COUNT=$6
NCPU=$7
echo "set vars"
echo $RAW_FILE
echo $TEMP_DIR
# create output directory, create file containing unannotated cdr3 sequences 
OUTPUT_LOCATION=$(Rscript scripts/igor_scripts/igor_preprocessing.R $RAW_FILE $TEMP_DIR)
echo "got location"
# submit job to run igor for specified individual, sample annotations, and reformat file
python scripts/igor_scripts/run_igor_custom.py $OUTPUT_LOCATION $OUTPUT_DIR $ANNOTATION_COUNT $NCPU $MODEL_PARAMS $MODEL_MARG
echo "finished"
# delete unnessary files
rm -r $OUTPUT_LOCATION
rm -r $RAW_FILE
echo "removed files"
