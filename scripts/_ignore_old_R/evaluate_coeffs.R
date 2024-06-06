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

MODEL_TYPE <<- args[7]

source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/model_fitting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/individual_comparison_functions.R'))

# Compile data for all subjects
motif_data = aggregate_all_subject_data(trim_type = TRIM_TYPE)

# bootstrap model fit
bootstrap_result = cluster_bootstrap_model_fit(motif_data, iter = 1000, trim_type = TRIM_TYPE, gene_type = GENE_NAME)
bootstrap_result$original_model_fit = FALSE

# fit model with original dataset
if (MODEL_TYPE %like% 'motif'){
    pwm_matrix = get_coefficient_matrix(motif_data, ref_base = 'A', trim_type = TRIM_TYPE)
    pwm_dt = as.data.table(pwm_matrix$result)
    setnames(pwm_dt, colnames(pwm_dt), map_chr(str_split(colnames(pwm_dt), '\\.', 2), 2))
    if (class(pwm_matrix) == 'list'){
        pwm_dt$base = rownames(pwm_matrix$result[[1]])
    } else {
        pwm_dt$base = rownames(pwm_matrix$result)
    }
    pwm_dt$snp_interaction = FALSE
    if (MODEL_TYPE %like% 'snp-interaction'){
        pwm_snp_dt = as.data.table(pwm_matrix$snp_interaction_result)
        setnames(pwm_snp_dt, colnames(pwm_snp_dt), map_chr(str_split(colnames(pwm_snp_dt), '\\.', 2), 2))
        if (class(pwm_matrix) == 'list'){
            pwm_snp_dt$base = rownames(pwm_matrix$snp_interaction_result[[1]])
        } else {
            pwm_snp_dt$base = rownames(pwm_matrix$snp_interaction_result)
        }
        pwm_snp_dt$snp_interaction = TRUE
        pwm_dt = rbind(pwm_dt, pwm_snp_dt)
    }
    model = pwm_matrix$model
} else {
    pwm_dt = NULL
    model = fit_model(motif_data, trim_type = TRIM_TYPE)
}

coefs = format_model_coefficient_output(model, pwm_dt, trim_type = TRIM_TYPE)
if (!('snp_interaction' %in% colnames(coefs))){
    coefs$snp_interaction = FALSE
}

coefs[snp_interaction == TRUE, parameter := paste0(parameter, ':snp')]
coefs$original_model_fit = TRUE

# write results
results = get_coef_pvalues(bootstrap_result, coefs)
filename =  get_model_bootstrap_file_name() 
fwrite(results, filename, sep = '\t')
