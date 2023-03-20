#!/bin/bash

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate mh_r 
set -eu

LOCUS=$1
TRIM_TYPE=$2
JOINING_GENE=$3
DATA_DIR=$4
NT_COUNT=$5
PRODUCTIVITY=$6

Rscript $PWD/analyses/different_processing_by_join.R $LOCUS $TRIM_TYPE $JOINING_GENE $DATA_DIR $NT_COUNT $PRODUCTIVITY 
