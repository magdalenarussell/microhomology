#!/bin/bash

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate mechanistic-trimming_py 
echo "set environ"
set -eu

RAW_FILE=$1
TEMP_DIR=$2
ADAPTIVE=$3
LOCUS=$4
OUTPUT_DIR=$5
MODEL_PARAMS=$6
MODEL_MARG=$7
ANNOTATION_COUNT=$8
NCPU=$9

echo "set vars"
echo $RAW_FILE
echo $TEMP_DIR

OUTPUT_LOCATION=$(Rscript scripts/igor_scripts/annotate/igor_preprocessing.R $RAW_FILE $TEMP_DIR $ADAPTIVE $LOCUS $NCPU)
echo "got location"

# submit job to run igor for specified individual, sample annotations, and reformat file
python scripts/igor_scripts/annotate/run_igor_annotation.py $OUTPUT_LOCATION $OUTPUT_DIR $ANNOTATION_COUNT $NCPU $MODEL_PARAMS $MODEL_MARG
echo "finished"

# delete unnessary files
rm -r $OUTPUT_LOCATION
rm -r $RAW_FILE
echo "removed files"
