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
library(vroom)
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

MODEL_TYPE <<- 'motif_two-side-base-count-beyond' 

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- 10

source(paste0(MOD_PROJECT_PATH,'scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'scripts/model_fitting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'plotting_scripts/plotting_functions.R'))

files = list.files(get_model_subsample_path(), full.names = TRUE, pattern = '*tsv') 
results = data.table(vroom(files))
results$log_coef = results$coefficient/log(10)
results[, q25 := quantile(log_coef, 0.25), by = .(parameter, base, prop)]
results[, q75 := quantile(log_coef, 0.75), by = .(parameter, base, prop)]
results[, median := quantile(log_coef, 0.5), by = .(parameter, base, prop)]
results[, variance := var(log_coef), by = .(parameter, base, prop)]

cols = c('parameter', 'base', 'prop', 'q25', 'q75', 'median', 'variance')
condensed = unique(results[, ..cols])

condensed$prop = log2(condensed$prop)
condensed$prop = factor(condensed$prop, levels = rev(unique(condensed$prop)))

names = data.table(parameter = c('left_base_count_GC', 'right_base_count_GC', 'right_base_count_AT', 'motif_5end_pos1', 'motif_3end_pos2', 'motif_3end_pos1'), nice_parameter = c('5\' GC base-count beyond', '3\' GC base-count beyond', '3\' AT base-count beyond', '5\' motif pos 1, ', '3\' motif pos 2, ', '3\' motif pos 1, '))

condensed = merge(condensed, names, by = 'parameter')
condensed[!is.na(base), nice_parameter := paste0(nice_parameter, base, ' nt')]
condensed$nice_parameter = factor(condensed$nice_parameter, levels = condensed[prop == -9][order(variance)]$nice_parameter)

label_data = condensed[prop == -9]
color_palette = set_color_palette(unique(condensed$nice_parameter))

require(ggrepel)
plot = ggplot(condensed[order(match(nice_parameter, levels(condensed$nice_parameter)))]) +
    geom_line(aes(x = prop, y = variance, color = nice_parameter, group = nice_parameter), size = 4, alpha = 0.8) +
    geom_point(aes(x = prop, y = variance, group = nice_parameter, color = nice_parameter), size = 8) +
    geom_text_repel(data = label_data, aes(y = variance, x = prop, label = nice_parameter, color = nice_parameter), nudge_x = 0.2, fontface = "bold", size = 11, direction = 'y', hjust = 0, point.padding = 1, max.overlaps = Inf, lineheight = 0.8) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('\nlog2(Proportion of full training data set)') +
    ylab('Variance of log10(probability of deletion)\n') +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    theme(legend.position = 'none', text = element_text(size = 45), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 40), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
    scale_x_discrete(expand = expansion(add = c(0.2, 4))) +
    scale_color_manual(values = color_palette)





path = get_manuscript_path()
file_name = paste0(path, '/subsample.pdf')
ggsave(file_name, plot = plot, width = 23, height = 18, units = 'in', dpi = 750, device = cairo_pdf)
