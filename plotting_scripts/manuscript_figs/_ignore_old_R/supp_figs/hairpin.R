source('config/config.R')

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

ANNOTATION_TYPE <<- 'igor' 

TRIM_TYPE <<- 'v_trim'

PRODUCTIVITY <<- 'nonproductive' 

MOTIF_TYPE <<- 'unbounded' 

NCPU <<- 2

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
stopifnot(GENE_NAME == 'v_gene')

MODEL_GROUP <<- 'all_subjects'

GENE_WEIGHT_TYPE <<- 'p_gene_marginal' 

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(1)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(2)

UPPER_TRIM_BOUND <<- as.numeric(14) 
LOWER_TRIM_BOUND <<- 2 

MODEL_TYPE <<- 'motif_two-side-base-count-beyond' 
stopifnot(MODEL_TYPE != 'null')

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- 10

TYPE <<- 'log_loss'
all_types = c('log_loss', 'expected_log_loss', 'v_gene_family_loss', 'log_loss_j_gene', 'full_v_gene_family_loss')

source(paste0(MOD_PROJECT_PATH,'plotting_scripts/plotting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'scripts/model_evaluation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'plotting_scripts/model_evaluation_functions.R'))

# get all model evaluation results (for all held-out groups)
all_eval_results = data.table()
for (type in all_types) {
    temp_eval_results = compile_evaluation_results(type)
    setnames(temp_eval_results, type, 'loss')
    temp_eval_results$loss_type = type
    all_eval_results = rbind(all_eval_results, temp_eval_results, fill = TRUE)
}

# specify held out cluster for "most different" sequence protocol
all_eval_results[loss_type %like% 'v_gene_family_loss', loss_type := paste0(loss_type, ', cluster ', held_out_clusters)]

# make names neat and assign color palette colors
motif_types = list.files(path = 'scripts/motif_class_functions/')
motif_types = str_sub(motif_types, end = -3)
nice_hairpin_names = make_hairpin_names_neat(motif_types) 
colors = set_color_palette(c(nice_hairpin_names))

#filter out extra "most different" held-out group trials
all_eval_results = all_eval_results[!(loss_type == "v_gene_family_loss, cluster 2" | loss_type == "full_v_gene_family_loss, cluster 2")]
all_eval_results = all_eval_results[!(loss_type == "v_gene_family_loss, cluster 3" | loss_type == "full_v_gene_family_loss, cluster 3")]
all_eval_results = all_eval_results[!(loss_type == "v_gene_family_loss, cluster 4" | loss_type == "full_v_gene_family_loss, cluster 4")]

# order loss types
loss_types = unique(all_eval_results$loss_type)
nice_loss_types = c('full V-gene\ntraining\ndataset', 'many held-out\nsubsets of\nV-gene\ntraining\ndataset', '\"most different\"\ncluster of\nV-genes\n(terminal seqs)', 'full J-gene\ndataset', '\"most different\"\ncluster of\nV-genes\n(full seqs)')
all_eval_results$loss_type = mapvalues(all_eval_results$loss_type, from = loss_types, to=nice_loss_types)

# make plot
plot = plot_model_evaluation_loss_paracoord(all_eval_results, model_type_list = MODEL_TYPE, left_motif_size_filter = LEFT_NUC_MOTIF_COUNT, right_motif_size_filter = RIGHT_NUC_MOTIF_COUNT, terminal_melting_5_end_length_filter = c(NA, LEFT_SIDE_TERMINAL_MELT_LENGTH), loss_bound = c(1.98, 2.26), color_palette = colors, same_motif_type = FALSE, custom_name = paste0(MODEL_TYPE, '_hairpin_nick_analysis'), plot_size = c(40, 25), write_plot = FALSE, expand_var = 2.2)

# save plot
path = get_manuscript_path()
file_name = paste0(path, '/hairpin_motif_base_count.pdf')
ggsave(file_name, plot = plot, width = 30, height = 20, units = 'in', dpi = 750, device = cairo_pdf)
