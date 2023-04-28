#!/bin/bash

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate mh_r 
set -eu

DATA_TYPE=$1
LOCUS=$2
TRIM_TYPE=$3
JOINING_GENE=$4
DATA_DIR=$5
NT_COUNT=$6
PRODUCTIVITY=$7
NCPU=$8
LOWER_BOUND=$9
UPPER_BOUND=${10}

Rscript $PWD/analyses/different_processing_by_join.R $DATA_TYPE $LOCUS $TRIM_TYPE $JOINING_GENE $DATA_DIR $NT_COUNT $PRODUCTIVITY $NCPU $LOWER_BOUND $UPPER_BOUND 
