#!/bin/bash

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate mechanistic-trimming_py 
echo "set environ"
set -eu

RAW_FILE_PATH=$1
TEMP_DIR=$2
ADAPTIVE=$3
SAMPLE_SIZE1=$4
SAMPLE_SIZE2=$5
LOCUS=$6
OUTPUT_DIR=$7
NCPU=$8

for RAW_FILE in $RAW_FILE_PATH/*.tsv; do
    OUTPUT_LOCATION=$(Rscript scripts/igor_scripts/igor_pretraining_preprocessing.R $RAW_FILE $TEMP_DIR $ADAPTIVE $SAMPLE_SIZE1 $LOCUS $NCPU)
    echo "finished processing for $RAW_FILE"
done

echo "got location"

# train igor model
python scripts/igor_scripts/train_new_model.py $OUTPUT_LOCATION $OUTPUT_DIR $SAMPLE_SIZE2
echo "finished"

# delete unnessary files
rm -r $OUTPUT_LOCATION
echo "removed files"
