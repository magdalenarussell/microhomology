source('config/config.R')
source('config/file_paths.R')

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(vroom)
library(Biostrings)
library(mclogit)
library(matrixcalc)
library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)

args = commandArgs(trailingOnly=TRUE)

PARAM_GROUP <<- 'nonproductive_v-j_trim_ligation-mh'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))
NCPU <<- as.numeric(args[1])
# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[2])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[3])
MODEL_TYPE <<- args[4]
L2 <<- args[5]
ANNOTATION_TYPE <<- args[6]

source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

# get full trained coefficients
coef_path = get_model_coef_file_path(L2 = L2)
full_coefs = fread(coef_path)

# get all trained coefficients
ex_boot_coef_path = get_bootstrap_model_coef_file_path(iteration = 1, L2 = L2)
boot_coef_paths = list.files(dirname(ex_boot_coef_path), full.names = TRUE, pattern = 'trained_coefs')
boot_coefs = data.table()
for (path in boot_coef_paths){
    temp = fread(path)
    boot_coefs = rbind(boot_coefs, temp)
}

# get bootstrap error
cols = c('coefficient', 'base', 'position', 'side', 'trim_type')
boot_coefs_cond = boot_coefs[, sd(value), by = cols]
setnames(boot_coefs_cond, 'V1', 'bootstrap_sd')

# combine data and compute zscores/pvalues
tog = merge(full_coefs, boot_coefs_cond)
tog[, zscore := value/bootstrap_sd]
tog[, pvalue := 2 * (1 - pnorm(abs(zscore)))]
test_count = nrow(tog)
tog[, adjusted_pvalue := pvalue * test_count]

# write results
path = get_bootstrap_pvalue_path(L2)
fwrite(tog[, -c('error')], path, sep = '\t')
