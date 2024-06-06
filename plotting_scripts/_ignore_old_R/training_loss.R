source('mechanistic-trimming/config/config.R')

library(cli)
library(devtools)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(foreach)
library(doParallel)
library(tidyverse)
library(plyr)
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

TRIM_TYPE <<- args[2]
trim_types = list.files(path = paste0(MOD_PROJECT_PATH, '/scripts/gene_specific_functions/'))
trim_types = str_sub(trim_types, end = -3)
stopifnot(TRIM_TYPE %in% trim_types)

PRODUCTIVITY <<- args[3]

MOTIF_TYPE <<- args[4] 
motif_types = list.files(path = paste0(MOD_PROJECT_PATH, '/scripts/motif_class_functions/'))
motif_types = str_sub(motif_types, end = -3)
stopifnot(MOTIF_TYPE %in% motif_types)

NCPU <<- as.numeric(args[5])

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')

MODEL_GROUP <<- 'all_subjects' 

GENE_WEIGHT_TYPE <<- args[6]
weight_types = list.files(path = paste0(MOD_PROJECT_PATH, '/scripts/sampling_procedure_functions/'))
weight_types = str_sub(weight_types, end = -3)
stopifnot(GENE_WEIGHT_TYPE %in% weight_types)

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[7])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[8])

UPPER_TRIM_BOUND <<- as.numeric(args[9]) 
LOWER_TRIM_BOUND <<- as.numeric(args[10])

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(args[11])

TYPE <<- 'log_loss'
all_types = c('log_loss', 'expected_log_loss', 'v_gene_family_loss', 'log_loss_j_gene', 'full_v_gene_family_loss')

LOSS_GENE_WEIGHT <<- 'p_gene_pooled'

source(paste0(MOD_PROJECT_PATH,'/scripts/model_evaluation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/model_evaluation_functions.R'))

# get all loss results
all_eval_results = data.table()
for (type in all_types) {
    temp_eval_results = compile_evaluation_results(type)
    setnames(temp_eval_results, type, 'loss')
    temp_eval_results$loss_type = type
    all_eval_results = rbind(all_eval_results, temp_eval_results, fill = TRUE)
}

# assign cluster for "most different" sequence protocol
all_eval_results[loss_type %like% 'v_gene_family_loss', loss_type := paste0(loss_type, ', cluster ', held_out_clusters)]

# get model types, and make them neat
maps = c("motif_two-side-base-count-beyond" = '1x2motif + two-side\nbase-count beyond\n(12 params)', 'null' = 'null (0 params)', 'motif' = "1x2motif (9 params)", 'two-side-base-count' = 'two-side base-count\n(3 params)', 'dna_shape-std' = '1x2DNA-shape\n(13 params)', 'linear-distance' = 'length (1 param)', 'motif_linear-distance' = '1x2motif + length\n(10 params)')

# pre-filter data
subset_eval_data = process_model_evaluation_file(all_eval_results, names(maps), left_motif_size_filter = LEFT_NUC_MOTIF_COUNT, right_motif_size_filter = RIGHT_NUC_MOTIF_COUNT, terminal_melting_5_end_length_filter = c(NA, LEFT_SIDE_TERMINAL_MELT_LENGTH))
    
maps = maps[names(maps) %in% unique(subset_eval_data$model_type)]

# reassign model names
subset_eval_data$model_type = mapvalues(subset_eval_data$model_type, from = names(maps), to = maps) 

# make model names neater and set palette colors
neat_names = make_model_names_neat(maps)
colors = set_color_palette(c(neat_names, '2x4motif (18 params)'), with_params = TRUE)

# get model evaluation results for the 2x4 motif model
eval_data_murugan = process_model_evaluation_file(all_eval_results, 'motif', 2, 4, NA)
eval_data_murugan$model_type = mapvalues(eval_data_murugan$model_type, from = 'motif', to = '2x4motif (18 params)')

# combine data sets
eval_tog = rbind(subset_eval_data, eval_data_murugan)

# filter out model evaluation types
eval_tog = eval_tog[!(loss_type %like% "v_gene_family_loss, cluster 1" | loss_type %like% "full_v_gene_family_loss, cluster 1")]
eval_tog = eval_tog[!(loss_type == "v_gene_family_loss, cluster 2" | loss_type == "full_v_gene_family_loss, cluster 2")]
eval_tog = eval_tog[!(loss_type == "v_gene_family_loss, cluster 4" | loss_type == "full_v_gene_family_loss, cluster 4")]

# order and reassign loss types (neater version)
loss_types = unique(eval_tog$loss_type)
# nice_loss_types = c('log_loss' = 'full V-gene\ntraining\ndataset', 'expected_log_loss' = 'many random,\nheld-out subsets\nof V-gene\ntraining\ndataset', '\"most different\"\ncluster of\nV-genes\n(terminal seqs)', 'full J-gene\ndataset', '\"most different\"\ncluster of\nV-genes\n(full seqs)')

loss_maps = c('log_loss' = 'full V-gene\ntraining\ndataset', 'expected_log_loss' = 'many random,\nheld-out subsets\nof V-gene\ntraining\ndataset', "v_gene_family_loss, cluster 1" = 'cluster 1 of V-genes\n(terminal seqs)', "v_gene_family_loss, cluster 3" = 'cluster 3 of V-genes\n(terminal seqs)', "v_gene_family_loss, cluster 1, 3" = 'cluster 1/3 of V-genes\n(terminal seqs)', "full_v_gene_family_loss, cluster 1" = 'cluster 1 of V-genes\n(full seqs)', "full_v_gene_family_loss, cluster 3" = 'cluster 3 of V-genes\n(full seqs)',"full_v_gene_family_loss, cluster 4" = 'cluster 4 of V-genes\n(full seqs)', "full_v_gene_family_loss, cluster 1, 3, 4" = 'cluster 1/3/4 of V-genes\n(full seqs)')

loss_maps = loss_maps[names(loss_maps) %in% loss_types]

eval_tog$loss_type = mapvalues(eval_tog$loss_type, from = names(loss_maps), to=loss_maps)

eval_tog = eval_tog[motif_type == MOTIF_TYPE]

# create plot
plot = plot_model_evaluation_loss_paracoord(eval_tog, model_type_list = c(new_model_types, '2x4motif'), pre_filter = TRUE, left_motif_size_filter = LEFT_NUC_MOTIF_COUNT, right_motif_size_filter = RIGHT_NUC_MOTIF_COUNT, terminal_melting_5_end_length_filter = c(NA, LEFT_SIDE_TERMINAL_MELT_LENGTH), loss_bound = c(1.92, 3.821), color_palette = colors, write_plot = FALSE, expand_var = 1.5) +
    ylab('Expected per-sequence log loss\n')

# save plot
path = get_manuscript_path()
file_name = paste0(path, '/loss_compare.pdf')
ggsave(file_name, plot = plot, width = 35, height = 18, units = 'in', dpi = 750, device = cairo_pdf)
