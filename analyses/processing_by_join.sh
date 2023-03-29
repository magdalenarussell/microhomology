#!/bin/bash

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate mh_r 
set -eu

LOCUS=$1
TRIM_TYPE=$2
JOINING_GENE=$3
DATA_DIR=$4
PRODUCTIVITY=$5

Rscript $PWD/analyses/sig_processing_by_join_multinom.R $LOCUS $TRIM_TYPE $JOINING_GENE $DATA_DIR $PRODUCTIVITY 
