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
LEFT_NUC_MOTIF_COUNT <<- as.numeric(1)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(2)

UPPER_TRIM_BOUND <<- as.numeric(14) 

MODEL_TYPE <<- 'motif_linear-distance_two-side-base-count-beyond-prop_snp-interaction-20717772'

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- 10

source(paste0(MOD_PROJECT_PATH,'scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'scripts/model_fitting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'plotting_scripts/plotting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'plotting_scripts/individual_comparison_functions.R'))

# get pwm from coefficient bootstrap
filename =  get_model_bootstrap_file_name() 
pwm = fread(filename)
pwm = pwm[parameter %like% 'snp']
pwm$parameter = str_replace(pwm$parameter, ':snp', '')
pwm$parameter = str_replace(pwm$parameter, '_prop', '')

# SNP genotype is opposite of convention, so switch sign of coefficient
pwm$coefficient = -1*(pwm$coefficient)

# plot coefficient heatmap
heatmap = plot_model_coefficient_heatmap_single_group(pwm, with_values = FALSE, write_plot = FALSE, limits = c(-0.007, 0.007))
heatmap = heatmap + 
    theme(legend.position = 'none', text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 20)) 

heatmap2 = plot_base_count_coefficient_heatmap_single_group(pwm, with_values = FALSE, write_plot = FALSE, limits = c(-0.007, 0.007))
heatmap2 = heatmap2 + 
    theme(legend.position = 'none', text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 20)) 

heatmap3 = plot_lindistance_coefficient_heatmap_single_group(pwm, with_values = FALSE, write_plot = FALSE, limits = c(-0.007, 0.007))
heatmap3 = heatmap3 + 
    theme(text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 18)) 

# align and combine plots
all = align_plots(heatmap, heatmap2, heatmap3, align = 'vh', axis = 'lb')

first_grid = plot_grid(all[[1]], all[[2]], nrow = 1, labels = c("A", "B"), label_size = 35, rel_widths = c(1, 1), align = 'h')
second_grid = plot_grid(all[[3]], NULL, ncol = 2, labels = c('C', ''), label_size = 35, align = 'h', rel_widths = c(1, 0.2))

together = plot_grid(first_grid, second_grid, nrow = 2, align = 'h')

# save plot
path = get_manuscript_path()
file_name = paste0(path, '/snp_interaction_heatmaps_motif_linear-distance_two-side-base-count-beyond.pdf')
ggsave(file_name, plot = together, width = 14, height = 12, units = 'in', dpi = 750, device = cairo_pdf)
