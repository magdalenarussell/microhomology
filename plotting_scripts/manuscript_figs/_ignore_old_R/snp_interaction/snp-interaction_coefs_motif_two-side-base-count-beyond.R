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

MODEL_TYPE <<- 'motif_two-side-base-count-beyond_snp-interaction-20717772'

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- 10

source(paste0(MOD_PROJECT_PATH, 'scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'scripts/model_fitting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'plotting_scripts/plotting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'plotting_scripts/individual_comparison_functions.R'))

# get coefficient bootsrap results
filename =  get_model_bootstrap_file_name() 
pwm = fread(filename)
pwm = pwm[parameter %like% 'snp']
pwm$parameter = str_replace(pwm$parameter, ':snp', '')
pwm$parameter = str_replace(pwm$parameter, '_prop', '')

# SNP genotype is opposite of convention, so switch sign of coefficient
pwm$coefficient = -1*(pwm$coefficient)

cols = colnames(pwm)[!(colnames(pwm) %in% c('snp_interaction', 'original_model_fit', 'zstat', 'iterations'))]
heatmap_data = pwm[, ..cols]
fwrite(heatmap_data, paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/snp_interaction/coefs.tsv'), sep = '\t')



# plot coefficient heatmaps
heatmap = plot_model_coefficient_heatmap_single_group(pwm, with_values = FALSE, write_plot = FALSE, limits = c(-0.011, 0.011))
heatmap = heatmap + 
    theme(legend.position = 'none', text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 20)) 

heatmap2 = plot_base_count_coefficient_heatmap_single_group(pwm, with_values = FALSE, write_plot = FALSE, limits = c(-0.011, 0.011))
heatmap2 = heatmap2 + 
    theme(text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 20),legend.position = 'bottom', legend.direction = 'horizontal', legend.justification="center") +
    guides(fill = guide_colourbar(barwidth = 35, barheight = 2)) +
    labs(fill = 'log10(probability of deletion)\t')

# extract legend
legend = get_legend(heatmap2) 
heatmap2 = heatmap2 + theme(legend.position = 'none')

# align and combine plots
all = align_plots(heatmap, heatmap2, legend, align = 'vh', axis = 'lbr')

first_grid = plot_grid(all[[1]], all[[2]], nrow = 1, labels = c("A", "B"), label_size = 35, rel_widths = c(1, 1), align = 'h')
first_grid = plot_grid(first_grid, NULL, all[[3]], nrow = 3, rel_heights = c(1, 0.02,0.2)) 

# save plot
path = get_manuscript_path()
file_name = paste0(path, '/snp_interaction_heatmaps.pdf')
ggsave(file_name, plot = first_grid, width = 18, height = 6, units = 'in', dpi = 750, device = cairo_pdf)


