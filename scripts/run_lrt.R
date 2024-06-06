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

ANNOTATION_TYPE <<- 'igor_sim_alpha'
PARAM_GROUP <<- 'nonproductive_v-j_trim_ligation-mh'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))
NCPU <<- as.numeric(args[1])
# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 1
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 2
MODEL_TYPE <<- 'motif_two-side-base-count-beyond_ligation-mh'
L2 <<- args[2]

source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

ANNOTATION_TYPE <<- args[3]

# Read processed data
filename = processed_data_path(sample_annotation = SAMPLE_ANNOT)
full_data = fread(filename)

full_pred_path = get_model_predictions_file_path(L2)
full_pred = fread(full_pred_path)


MODEL_TYPE2 <<- 'motif_two-side-base-count-beyond'
sub_pred_path = get_model_predictions_file_path(L2, model_type=MODEL_TYPE2)
sub_pred = fread(sub_pred_path)

# get log likelihood
full_pred[, seq_log_lik := count * log(predicted_prob)]
full_pred[, logLik := sum(seq_log_lik)]

sub_pred[, seq_log_lik := count * log(predicted_prob)]
sub_pred[, logLik := sum(seq_log_lik)]

# convert loss to likelihood
total_seq_count = sum(full_data$count)

# LRT test statistic
D = -2 * (unique(sub_pred$logLik) - unique(full_pred$logLik))

# LRT pval 
p_value = pchisq(D, 1, lower.tail = FALSE)

# compile results
df = data.table(annotation_type = ANNOTATION_TYPE, 
                full_model = MODEL_TYPE, 
                sub_model = MODEL_TYPE2,
                full_model_logLik = unique(full_pred$logLik),
                sub_model_logLik = unique(sub_pred$logLik),
                LRT_test_statistic = D,
                LRT_pval = p_value)

# write results
lrt_path = get_LRT_file_path(L2)

if (file.exists(lrt_path)){
    lrt_df = fread(lrt_path)
    lrt_df = rbind(lrt_df, df)
    fwrite(lrt_df, lrt_path, sep = '\t')
} else {
    fwrite(df, lrt_path, sep = '\t')
}
