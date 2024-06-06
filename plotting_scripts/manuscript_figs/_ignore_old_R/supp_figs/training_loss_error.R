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

ANNOTATION_TYPE <<- 'igor' 

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

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(1)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(2)

UPPER_TRIM_BOUND <<- as.numeric(14) 
LOWER_TRIM_BOUND <<- 2 

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(10)

TYPE <<- 'expected_log_loss'

LOSS_GENE_WEIGHT <<- 'p_gene_given_subject'

source(paste0(MOD_PROJECT_PATH,'scripts/model_evaluation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'plotting_scripts/plotting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'plotting_scripts/model_evaluation_functions.R'))

# get all loss results
path_name = get_per_run_model_evaluation_path(TYPE)
path_name = str_replace(path_name, '/expected_log_loss', '/expected_log_loss/raw')

require(vroom)
files = list.files(path_name, full.names = TRUE)
all = vroom(files)
all = as.data.table(all)

motif24 = all[model_type == 'motif' & motif_length_5_end == 2]
all = all[!(model_type == 'motif' & motif_length_5_end == 2)]

# get model types, and make them neat
orig_model_types = c("motif_two-side-base-count-beyond", "null", "motif", "two-side-base-count", 'dna_shape-std', 'linear-distance', 'motif_linear-distance')
new_model_types = c('1x2motif + two-side\nbase-count beyond\n(12 params)', 'null (0 params)', "1x2motif (9 params)", 'two-side base-count\n(3 params)', '1x2DNA-shape\n(13 params)', 'length (1 param)', '1x2motif + length\n(10 params)')

# reassign model names
all$model_type = mapvalues(all$model_type, from = orig_model_types, to = new_model_types) 

# make model names neater and set palette colors
neat_names = make_model_names_neat(new_model_types)
colors = set_color_palette(c(neat_names, '2x4motif (18 params)'), with_params = TRUE)

motif24$model_type = '2x4motif (18 params)'

all = rbind(all, motif24)
all[, mean := mean(expected_log_loss), by = .(motif_length_5_end, model_type)]

order = unique(all[order(mean)]$model_type)
all$model_type = factor(all$model_type, levels = order)

plot = ggplot(all) +
    geom_density(aes(x = expected_log_loss, color = model_type, fill = model_type), alpha = 0.6) +
    geom_vline(aes(xintercept = mean, color = model_type), size = 5) +
    facet_grid(rows = vars(model_type)) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('\nExpected per-sequence log loss') +
    ylab('Density\n')+
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    theme(legend.position = 'none', text = element_text(size = 44), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 44), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), strip.text.y = element_text(angle = 0)) + 
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors)



# save plot
path = get_manuscript_path()
file_name = paste0(path, '/expected_loss_error.pdf')
ggsave(file_name, plot = plot, width = 32, height = 24, units = 'in', dpi = 750, device = cairo_pdf)
