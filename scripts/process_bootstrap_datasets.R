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

# Read processed data
filename = processed_data_path(sample_annotation = SAMPLE_ANNOT)
full_df = fread(filename)

# bootstrap data
set.seed(1)
total = 100
subset = unique(full_df[!is.na(count)][, c('index', 'count')])

for (i in seq(total)){
    obs_count = sum(subset$count, na.rm = TRUE)
    sampled_indices = sample(x = subset$index, size = obs_count, replace = TRUE, prob = subset$count)
    sample_df = data.table(index = sampled_indices)
    sample_df = sample_df[, .N, by = index]
    setnames(sample_df, 'N', 'count')

    sample_tog = merge(full_df[, -c('count', 'weighted_observation', 'total_tcr')], sample_df, by = c('index'), all = TRUE)
    sample_tog[is.na(count) & (frame_stop == TRUE | frame_type == 'Out'), count := 0]
    sample_tog[is.na(count) & (frame_stop == FALSE & frame_type == 'In'), count := NA]

    sample_tog_proc = calculate_subject_gene_weight(sample_tog, gene_type = GENE_NAME, trim_type = TRIM_TYPE, only_nonprod_sites = ONLY_NONPROD_SITES, sample_annotation = SAMPLE_ANNOT)
    sample_tog_proc = subset_processed_data(sample_tog_proc, trim_type = TRIM_TYPE, gene_type = GENE_NAME)

    boot_name = bootstrap_data_path(iteration = i, sample_annotation=SAMPLE_ANNOT)
    fwrite(sample_tog_proc, boot_name, sep = '\t')
}
