#!/bin/bash

set -eu

OUTPUT_DIR=$1
INPUT_DIR=$2
SEQ_COUNT=$3
NCPU=$4
PARTITION=$5
CHAIN=$6

for ITER in {1..1000..1}; do
    COMMAND="sbatch -c $NCPU -p $PARTITION -q $PARTITION scripts/simulate_seqs_with_igor.sh $OUTPUT_DIR $INPUT_DIR $SEQ_COUNT $NCPU $CHAIN $ITER"
    echo "Running $COMMAND"
    $COMMAND
done
