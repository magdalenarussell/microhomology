source('mechanistic-trimming/config/config.R')

library(ggplot2)
library(cowplot)
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


ANNOTATION_TYPE<<- 'igor' 
TRIM_TYPE <<- 'v_trim' 
trim_types = list.files(path = 'scripts/gene_specific_functions/')
trim_types = str_sub(trim_types, end = -3)
stopifnot(TRIM_TYPE %in% trim_types)

PRODUCTIVITY <<- 'nonproductive' 

MOTIF_TYPE <<- 'unbounded' 
motif_types = list.files(path = 'scripts/motif_class_functions/')
motif_types = str_sub(motif_types, end = -3)
stopifnot(MOTIF_TYPE %in% motif_types)

NCPU <<- as.numeric(2)

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')

MODEL_GROUP <<- 'all_subjects'

GENE_WEIGHT_TYPE <<- 'p_gene_marginal' 
stopifnot(GENE_WEIGHT_TYPE %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform'))

LOWER_TRIM_BOUND <<- 2
UPPER_TRIM_BOUND <<- 14 

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- 10 

TYPE <<- 'log_loss' 
stopifnot(TYPE %in% c('log_loss', 'expected_log_loss', 'aic', 'raw_loss', 'old_loss_cv', 'log_loss_j_gene', 'v_gene_family_loss', 'full_v_gene_family_loss'))


source(paste0(MOD_PROJECT_PATH,'scripts/model_evaluation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'plotting_scripts/plotting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'plotting_scripts/model_evaluation_functions.R'))

# compile loss results for training data set
eval_results = compile_evaluation_results(TYPE)

# assign result for two-side-base-count to be 0-0 motif-base-count model 
model_type = 'motif_two-side-base-count-beyond' 
base_model = 'two-side-base-count'
eval_results[model_type == base_model & motif_type == MOTIF_TYPE, model_type := 'motif_two-side-base-count-beyond']

# plot loss heatmap
plot = plot_model_evaluation_heatmap(eval_results, TYPE, model_type = model_type, terminal_melting_5_end_length_filter= 10, write_plot = FALSE)

# get loss results for J-genes from training data set
TYPE <<- 'log_loss_j_gene' 

source('scripts/model_evaluation_functions.R')
source('plotting_scripts/plotting_functions.R')
source('plotting_scripts/model_evaluation_functions.R')

# compile J-gene loss results
eval_results2 = compile_evaluation_results(TYPE)

# assign result for two-side-base-count to be 0-0 motif-base-count model 
eval_results2[model_type == base_model & motif_type == MOTIF_TYPE, model_type := 'motif_two-side-base-count-beyond']

# plot loss heatmap
plot2 = plot_model_evaluation_heatmap(eval_results2, TYPE, model_type = model_type, terminal_melting_5_end_length_filter= 10, write_plot = FALSE)

# combine heatmaps
final = plot_grid(plot, NULL, plot2, ncol = 3, labels = c("A", "", "B"), label_size = 35, rel_widths = c(1, 0.05, 0.96), align = 'h')

# save plot
path = get_manuscript_path()
file_name = paste0(path, '/motif_base_count_size.pdf')
ggsave(file_name, plot = final, width = 20, height = 7, units = 'in', dpi = 750, device = cairo_pdf)
