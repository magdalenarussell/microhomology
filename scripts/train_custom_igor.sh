#!/bin/bash

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate mechanistic-trimming_py 
echo "set environ"
set -eu

RAW_FILE_PATH=$1
TEMP_DIR=$2
ADAPTIVE=$3
LOCUS=$4
SAMPLE_SIZE=$5
LOCUS=$6
OUTPUT_DIR=$7
NCPU=$8

echo "set vars"
echo $TEMP_DIR

for RAW_FILE in $RAW_FILE_PATH/*.tsv; do
    Rscript scripts/igor_scripts/train/igor_pretraining_preprocessing.R $RAW_FILE $TEMP_DIR $ADAPTIVE $SAMPLE_SIZE $LOCUS
done

python scripts/igor_scripts/train/train_new_model.py $TEMP_DIR $OUTPUT_DIR $SAMPLE_SIZE $LOCUS
echo "finished"

# delete unnessary files
rm -r $TEMP_DIR
echo "removed files"
