source('config/config.R')
source('config/file_paths.R')

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

PARAM_GROUP <<- args[1]
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))
NCPU <<- as.numeric(args[2])
# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[3])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[4])
MODEL_TYPE <<- args[5]
MODEL_TYPE2 <<- args[6]
L2 <<- args[7]
ANNOTATION_TYPE <<- args[8]
DF <<- as.numeric(args[9])

source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

# Read processed data
filename = processed_data_path(sample_annotation = SAMPLE_ANNOT)
full_data = fread(filename)

full_pred_path = get_model_predictions_file_path(L2)

full_pred = fread(full_pred_path)
if ('new_type' %in% colnames(full_pred)){
    full_pred = full_pred[new_type != '0']
}

sub_pred_path = get_model_predictions_file_path(L2, model_type=MODEL_TYPE2)
sub_pred = fread(sub_pred_path)
if ('new_type' %in% colnames(sub_pred)){
    sub_pred = sub_pred[new_type != '0']
}

if (SAMPLE_ANNOT == FALSE){
    full_pred[, index_prob_sum := sum(predicted_prob, na.rm = TRUE), by = index]
    full_pred[, count_prob :=  predicted_prob/index_prob_sum]
    full_pred[, adjusted_count := count * count_prob]

    sub_pred[, index_prob_sum := sum(predicted_prob, na.rm = TRUE), by = index]
    sub_pred[, count_prob :=  predicted_prob/index_prob_sum]
    sub_pred[, adjusted_count := count * count_prob]
} else {
    full_pred[, adjusted_count := count]
    sub_pred[, adjusted_count := count]
}

# get log likelihood
full_pred[, seq_log_lik := adjusted_count * log(predicted_prob)]
full_pred[, logLik := sum(seq_log_lik, na.rm = TRUE)]

sub_pred[, seq_log_lik := adjusted_count * log(predicted_prob)]
sub_pred[, logLik := sum(seq_log_lik, na.rm = TRUE)]

# LRT test statistic
D = -2 * (unique(sub_pred$logLik) - unique(full_pred$logLik))

# LRT pval 
p_value = pchisq(D, DF, lower.tail = FALSE)

# compile results
df = data.table(annotation_type = ANNOTATION_TYPE, 
                full_model = MODEL_TYPE, 
                sub_model = MODEL_TYPE2,
                full_model_logLik = unique(full_pred$logLik),
                sub_model_logLik = unique(sub_pred$logLik),
                LRT_test_statistic = D,
                LRT_pval = p_value,
                DF = DF)

# write results
lrt_path = get_LRT_file_path(L2, sample_annotation = SAMPLE_ANNOT)

if (file.exists(lrt_path)){
    lrt_df = fread(lrt_path)
    lrt_df = rbind(lrt_df, df)
    fwrite(lrt_df, lrt_path, sep = '\t')
} else {
    fwrite(df, lrt_path, sep = '\t')
}
