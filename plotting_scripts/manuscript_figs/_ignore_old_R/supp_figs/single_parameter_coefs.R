source('config/config.R')

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(Biostrings)
library(ggplot2)
library(cowplot)
library(mclogit)
library(matrixcalc)
library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)

args = commandArgs(trailingOnly=TRUE)

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

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(0)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(0)

UPPER_TRIM_BOUND <<- as.numeric(14) 

MODEL_TYPE <<- 'two-side-base-count' 

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- 10 

source(paste0(MOD_PROJECT_PATH,'scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'scripts/model_fitting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'plotting_scripts/plotting_functions.R'))

# Read in model coefficient data 
pwm = get_model_coefficient_data() 

# plot coefficient heatmap
heatmap = plot_base_count_coefficient_heatmap_single_group(pwm, with_values = FALSE, write_plot = FALSE)
heatmap = heatmap + 
    theme(text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 25)) 

# save plot
path = get_manuscript_path()
file_name = paste0(path, '/single_param_coefs.pdf')
ggsave(file_name, plot = heatmap, width = 16, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
