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

PARAM_GROUP <<- 'nonproductive_v-j_trim_ligation-mh'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))

NCPU <<- as.numeric(2)

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(1)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(2)

MODEL_TYPE <<- 'motif_two-side-base-count-beyond_average-mh_ligation-mh'

L2 <<- 'True'

ANNOTATION_TYPE <<- 'igor_alpha'

source(paste0(MOD_PROJECT_PATH, '/config/file_paths.R'))
source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

file_root = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/motif_base_count_average-mh_ligation-mh_model/', PARAM_GROUP)

# Read in model coefficient data 
coef_path = get_model_coef_file_path(L2)
coefs = fread(coef_path)

fwrite(coefs, paste0(file_root, '/coefs.tsv'), sep = '\t')

v_motif_heatmap = plot_motif_coefficient_heatmap_single_group(coefs[trim_type == 'v_trim'], with_values = FALSE, limits = c(-0.450, 0.450)) + ggtitle('   V-trimming coefficients')
j_motif_heatmap = plot_motif_coefficient_heatmap_single_group(coefs[trim_type == 'j_trim'], with_values = FALSE, limits = c(-0.450, 0.450)) + ggtitle('   J-trimming coefficients')

v_motif_heatmap = v_motif_heatmap + theme(legend.position = 'none') 
j_motif_heatmap = j_motif_heatmap + theme(legend.position = 'none') 

# plot base count heatmap of coefficients
v_base_count_heatmap = plot_base_count_coefficient_heatmap_single_group(coefs[trim_type == 'v_trim'], with_values = FALSE, limits = c(-0.450, 0.450))
j_base_count_heatmap = plot_base_count_coefficient_heatmap_single_group(coefs[trim_type == 'j_trim'], with_values = FALSE, limits = c(-0.450, 0.450))

v_base_count_heatmap = v_base_count_heatmap + theme(legend.position = 'none') 
j_base_count_heatmap = j_base_count_heatmap + theme(legend.position = 'none') 

# plot mh heatmap
mh_heatmap = plot_ligation_mh_coefficient_heatmap_single_group(coefs, with_values = FALSE, limits = c(-0.450, 0.450))
config_heatmap = plot_average_mh_coefficient_heatmap_single_group(coefs, with_values = FALSE, limits = c(-0.450, 0.450))

config_heatmap = config_heatmap + theme(legend.position = 'none')

# isolate legend
legend = get_legend(mh_heatmap) 
mh_heatmap = mh_heatmap + theme(legend.position = 'none')

all = align_plots(v_motif_heatmap, j_motif_heatmap, v_base_count_heatmap, j_base_count_heatmap, config_heatmap, mh_heatmap, legend, align = 'vh', axis = 'lbr')

first_grid = plot_grid(all[[1]], all[[2]], nrow = 1, rel_widths = c(1, 1), align = 'h')
second_grid = plot_grid(all[[3]], all[[4]], nrow = 1, rel_widths = c(1, 1), align = 'h')

third_grid = plot_grid(NULL, all[[5]], NULL, nrow = 1, rel_widths = c(0.1, 0.4, 0.1), align = 'h')
fourth_grid = plot_grid(NULL, all[[6]], NULL, nrow = 1, rel_widths = c(0.1, 0.4, 0.1), align = 'h')

assign('together', plot_grid(first_grid, NULL, second_grid, NULL, third_grid, NULL, fourth_grid, NULL, legend, nrow = 9, rel_heights = c(1, 0.08, 1, 0.08, 0.75, 0.08, 0.75, 0.08, 0.2)))

# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/motif_base_count_average-mh_ligation-mh_model/', PARAM_GROUP, '/coef_heatmap.pdf')

ggsave(file_name, plot = together, width = 18, height = 24, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)
