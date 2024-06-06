#!/bin/bash

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate mechanistic-trimming_py

set -eu

module load Automake/1.16.1-GCCcore-7.3.0
module load GCCcore/7.3.0

OUTPUT_DIR=$1
INPUT_DIR=$2
SEQ_COUNT=$3
NCPU=$4
CHAIN=$5
ITER=$6


cd $HOME/bin

TEMP_DIR=$INPUT_DIR/temp_${ITER}
mkdir $TEMP_DIR
#Now generate synthetic sequences from the provided human beta chain model
#This will create the directory bar_generate with the corresponding files containing the generated sequences and their realizations
igor -threads $NCPU -set_wd $TEMP_DIR -batch bar -species human -chain $CHAIN -generate $SEQ_COUNT 

# now annotate sequences 
cd $HOME/microhomology

COMMAND="python scripts/igor_scripts/simulate/convert_seqs_with_igor_${CHAIN}.py $OUTPUT_DIR $TEMP_DIR $ITER"

$COMMAND
