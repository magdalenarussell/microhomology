#!/bin/bash

#SBATCH --job-name=igor
#SBATCH -p compute 
#SBATCH -A stf
#SBATCH --nodes=1
#SBATCH --mem=170GB
#SBATCH --ntasks-per-node=40
#SBATCH --time=24:00:00

source ~/.bashrc

WDPATH=$1
INFILE=$2
BATCHNAME=$3
echo "WDPATH ${WDPATH}"
echo "INFILE ${INFILE}"
echo "BATCHNAME ${BATCHNAME}"

IGOR="igor -set_wd ${WDPATH} -batch ${BATCHNAME} -species human -chain alpha"

${IGOR} -read_seqs ${INFILE} #Read seqs
${IGOR} -align --all #  Align
${IGOR} -infer --L_thresh "1e-300" --N_iter 10 #  Infer
