#!/bin/bash

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate mechanistic-trimming_py 

set -eu

RAW_FILE=$1
TEMP_DIR=$2
OUTPUT_DIR=$3
ANNOTATION_COUNT=$4
NCPU=$5

# create output directory, create file containing unannotated cdr3 sequences 
OUTPUT_LOCATION=$(Rscript scripts/igor_scripts/igor_preprocessing.R $RAW_FILE $TEMP_DIR)

# submit job to run igor for specified individual
bash scripts/igor_scripts/run_igor.sh $OUTPUT_LOCATION $ANNOTATION_COUNT $NCPU

# sample from possible annotations
Rscript scripts/igor_scripts/igor_sample_annotations.R $OUTPUT_LOCATION $NCPU

# translate igor outputs to human readable
python scripts/igor_scripts/igor_output_processing.py $OUTPUT_LOCATION $OUTPUT_LOCATION

# reformat file
Rscript scripts/igor_scripts/igor_reformat.R $OUTPUT_LOCATION $OUTPUT_DIR

# delete unnessary files
rm -r $OUTPUT_LOCATION
