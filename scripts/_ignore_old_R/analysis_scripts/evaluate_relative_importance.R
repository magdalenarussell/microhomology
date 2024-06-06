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
stopifnot(ANNOTATION_TYPE %in% c('igor', 'parsimony', 'alpha'))

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

MODEL_TYPE <<- 'motif_two-side-base-count-beyond'

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(args[7])

if (!(grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE))){
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA
}

VALIDATION_DATA_DIR <<- args[8]
VALIDATION_TYPE <<- args[9]
VALIDATION_TRIM_TYPE <<- args[10]
VALIDATION_PRODUCTIVITY <<- args[11]
VALIDATION_GENE_NAME <<- paste0(substring(VALIDATION_TRIM_TYPE, 1, 1), '_gene')
stopifnot(VALIDATION_TYPE %in% c('validation_data_alpha', 'validation_data_beta', 'validation_data_gamma', 'validation_data_delta', 'validation_data_igh', 'validation_data_igk', 'validation_data_igl', 'igor'))

LOSS_GENE_WEIGHT <<- 'p_gene_pooled' 

source(paste0(MOD_PROJECT_PATH,'/scripts/model_fitting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/analysis_scripts/rel_importance_functions.R'))

model = load_model()
pwm = get_model_coefficient_data()
output_file_name = get_per_run_rel_importance_model_file_name()

ANNOTATION_TYPE <<- VALIDATION_TYPE
TYPE <<- 'validation_data' 
TRIM_TYPE <<- VALIDATION_TRIM_TYPE
GENE_NAME <<- VALIDATION_GENE_NAME
PRODUCTIVITY <<- VALIDATION_PRODUCTIVITY

source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH, '/scripts/model_fitting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/model_evaluation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/analysis_scripts/pwm_profile_functions.R'))

validation_data = aggregate_validation_data(directory = VALIDATION_DATA_DIR, trim_type = TRIM_TYPE, gene_type = GENE_NAME)
validation_data = get_model_feature_scores(validation_data, pwm, trim_type = TRIM_TYPE, gene_type = GENE_NAME)

rel_importance_model = fit_rel_importance_model(validation_data, gene_type = GENE_NAME)
loss = evaluate_loss(validation_data, rel_importance_model, trim_type = TRIM_TYPE, gene_type = GENE_NAME)
write_rel_importance_result_dt(coef(rel_importance_model), loss, output_file_name)
