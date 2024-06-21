#!/bin/bash

set -eu

module load Automake/1.16.1-GCCcore-7.3.0
module load GCCcore/7.3.0

OUTPUT_DIR=$1
ANNOTATION_COUNT=$2
NCPU=$3
CHAIN=$4

cd $HOME/bin

#We first read the sequences contained in a text file inside the demo folder
#This will create the align folder in the working directory and the mydemo_indexed_seqs.csv file.
igor -threads $NCPU -set_wd $OUTPUT_DIR -batch dataIGoRD4KcOPYTXn -read_seqs $OUTPUT_DIR/cdr3_seqs.txt 

#Now let's align the sequences against the provided human beta chain genomic templates with default parameters
#This will create foo_V_alignments.csv, foo_D_alignments.csv and foo_J_alignments.csv files inside the align folder.
igor -threads $NCPU -set_wd $OUTPUT_DIR -batch dataIGoRD4KcOPYTXn -species human -chain $CHAIN -align --all

# this will train a model
igor -threads $NCPU -set_wd $OUTPUT_DIR -batch dataIGoRD4KcOPYTXn -species human -chain $CHAIN -infer 
