source('config/config.R')

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(Biostrings)
library(mclogit)
library(matrixcalc)
library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)

args = commandArgs(trailingOnly=TRUE)

ANNOTATION_TYPE <<- args[1]
annotation_types = list.files(path = paste0(MOD_PROJECT_PATH, '/scripts/annotation_specific_functions/'))
annotation_types = str_sub(annotation_types, end = -3)
stopifnot(ANNOTATION_TYPE %in% annotation_types)

PARAM_GROUP <<- args[2]
param_types = list.files(path = paste0(MOD_PROJECT_PATH, '/scripts/param_groups/'))
param_types = str_sub(param_types, end = -3)
stopifnot(PARAM_GROUP %in% param_types)
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))

NCPU <<- as.numeric(args[3])

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[4])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[5])

MODEL_TYPE <<- args[6]

source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

# Compile data for all subjects
motif_data = aggregate_all_subject_data(trim_type = TRIM_TYPE, sample_annotation = SAMPLE_ANNOT, only_nonprod_sites=ONLY_NONPROD_SITES)

# subset data
motif_data = subset_processed_data(motif_data, trim_type = TRIM_TYPE, gene_type = GENE_NAME)

# Write processed data
filename = processed_data_path(sample_annotation = SAMPLE_ANNOT)
fwrite(motif_data, filename, sep = '\t')
