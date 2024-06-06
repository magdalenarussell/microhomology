source('mechanistic-trimming/config/config.R')

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
stopifnot(ANNOTATION_TYPE %in% c('igor_alpha', 'igor_beta', 'adaptive_alpha', 'adaptive_beta'))

PARAM_GROUP <<- args[2]
param_types = list.files(path = paste0(MOD_PROJECT_PATH, '/scripts/param_groups/'))
param_types = str_sub(param_types, end = -3)
stopifnot(PARAM_GROUP %in% param_types)
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))

NCPU <<- as.numeric(args[3])

# NOTE: This method is only applicable for models fit across all subjects!
MODEL_GROUP <<- 'all_subjects'

GENE_WEIGHT_TYPE <<- args[4]
weight_types = list.files(path = paste0(MOD_PROJECT_PATH, '/scripts/sampling_procedure_functions/'))
weight_types = str_sub(weight_types, end = -3)
stopifnot(GENE_WEIGHT_TYPE %in% weight_types)

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[5])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[6])

MODEL_TYPE <<- args[7]

source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/model_fitting_functions.R'))

# Compile data for all subjects
motif_data = aggregate_all_subject_data(trim_type = TRIM_TYPE)

# Fit model, write predicted distribution and pwm files
fit_model_by_group(motif_data, trim_type = TRIM_TYPE, gene_type = GENE_NAME)
model = fit_model(motif_data, formula = get_model_formula(trim_type = TRIM_TYPE, gene_type = GENE_NAME), trim_type = TRIM_TYPE)
save_model(model)
