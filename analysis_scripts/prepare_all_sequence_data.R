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

ANNOTATION_TYPE <<- 'igor_alpha'
PARAM_GROUP <<- args[1]
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))
NCPU <<- as.numeric(args[2])
# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 1
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 2
MODEL_TYPE <<- 'motif_two-side-base-count-beyond_average-mh_ligation-mh'

source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

# get IGoR Vtrim and Jtrim probs
jt = fread(paste0(MOD_OUTPUT_PATH, '/igor_alpha/nonproductive_v-j_trim_ligation-mh/unbounded_motif_trims_bounded_-2_14/1_2_motif_two-side-base-count-beyond_average-mh_ligation-mh/igor_prob_experiment/igor_alpha_j_trim_params.tsv'))
vt = fread(paste0(MOD_OUTPUT_PATH, '/igor_alpha/nonproductive_v-j_trim_ligation-mh/unbounded_motif_trims_bounded_-2_14/1_2_motif_two-side-base-count-beyond_average-mh_ligation-mh/igor_prob_experiment/igor_alpha_v_trim_params.tsv'))
jt$indicator = 1
vt$indicator = 1

tog = merge(vt, jt, by = "indicator", allow.cartesian = TRUE)
tog[, igor_joint_prob := v_trim_prob * j_trim_prob]
tog[, vj_insert := 0]

# prepare data for model predictions
processed = compile_data_for_subject(file_path=NULL, dataset=tog, write = FALSE)

# final model preparation, getting model parameter values for each obs
weighted_together = inner_aggregation_processing(processed, trim_type = TRIM_TYPE, gene_type = GENE_NAME)

# subset columns 
trim_vars = get_trim_order(TRIM_TYPE)
genes = get_gene_order(GENE_NAME)
params = get_parameter_vector(trim_vars, genes)
cols = c('index', paste0(genes), trim_vars, 'weighted_observation', 'count', 'total_tcr', params, 'frame_type', 'frame_stop')

df = weighted_together[, ..cols]

# re-normalize igor probabilities
cols2 = c(genes)
tog[, igor_prob_total := sum(igor_joint_prob), by = cols2]
tog[, norm_default_igor_prob := igor_joint_prob/igor_prob_total]

# combine igor default probs
igor_cols = c(genes, trim_vars, 'igor_joint_prob', 'norm_default_igor_prob')
df = merge(df, unique(tog[, ..igor_cols]), by = c(genes, trim_vars))

# get path
processed_path = processed_data_path()
path = file.path(dirname(processed_path), 'igor_prob_experiment')
dir.create(path, recursive = TRUE, showWarnings = FALSE)
filename = file.path(path, basename(processed_path))

fwrite(df, filename, sep = '\t')
cat(filename)
