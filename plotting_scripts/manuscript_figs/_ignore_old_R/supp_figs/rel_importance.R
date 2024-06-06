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
GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
PRODUCTIVITY <<- 'nonproductive'

MOTIF_TYPE <<- 'unbounded' 
motif_types = list.files(path = 'scripts/motif_class_functions/')
motif_types = str_sub(motif_types, end = -3)
stopifnot(MOTIF_TYPE %in% motif_types)

NCPU <<- as.numeric(2)

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

TYPE <<- 'log_loss'
LOSS_GENE_WEIGHT <<- 'p_gene_given_subject' 
stopifnot(LOSS_GENE_WEIGHT %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform', 'p_gene_marginal_all_seqs'))

source(paste0(MOD_PROJECT_PATH,'scripts/model_evaluation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'plotting_scripts/plotting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'plotting_scripts/model_evaluation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'scripts/analysis_scripts/rel_importance_functions.R'))

all_eval_results = data.table()
rel_import_results = data.table()

# set all validation, trim-type, and productivity variables
annotation_types = c('validation_data_alpha', 'validation_data_beta', 'igor', 'validation_data_gamma', 'validation_data_igh')
trim_types = c('j_trim', 'v_trim')
prods = c('productive', 'nonproductive')

# get all rel importance results
for (trim_type in trim_types){
    for (type in annotation_types) {
        for (prod in prods){
            VALIDATION_TYPE <<- type
            VALIDATION_PRODUCTIVITY <<- prod
            VALIDATION_TRIM_TYPE <<- trim_type

            source('scripts/analysis_scripts/rel_importance_functions.R')

            temp_rel_results = compile_rel_importance_results()
            temp_rel_results$loss_type = type
            temp_rel_results$trim_type = trim_type
            temp_rel_results$productivity = prod
            rel_import_results = rbind(rel_import_results, temp_rel_results, fill = TRUE)
        }
    }
}
setnames(rel_import_results, 'loss', 'new_model_loss')

# get all validation loss results
for (trim_type in trim_types){
    for (type in annotation_types) {
        for (prod in prods){
            
            ANNOTATION_TYPE <<- type
            TRIM_TYPE <<- trim_type
            GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
            PRODUCTIVITY <<- prod

            source('scripts/model_evaluation_functions.R')
            source('plotting_scripts/plotting_functions.R')
            source('plotting_scripts/model_evaluation_functions.R')
            temp_eval_results = compile_evaluation_results(type)
            setnames(temp_eval_results, type, 'loss')
            temp_eval_results$loss_type = type
            temp_eval_results$trim_type = trim_type
            temp_eval_results$productivity = prod
            all_eval_results = rbind(all_eval_results, temp_eval_results, fill = TRUE)
        }
    }
}

# combine data sets
together = merge(all_eval_results[model_type == 'motif_two-side-base-count-beyond'], rel_import_results, by = c('loss_type', 'trim_type', 'productivity'))

# convert loss names to long, neat version
together[, long_loss_type := paste0(loss_type, '_', productivity, '_', trim_type)]
loss_types = unique(together$long_loss_type)
nice_loss_types = c(rep('TCRB training dataset genes\n(not included during model training)', 2), 'TCRB training dataset', 'TCRB training dataset genes\n(not included during model training)', rep('TCRA testing dataset', 4), rep('TCRB testing dataset', 4), rep('TCRG testing dataset', 4), rep('IGH testing dataset', 4))
loss_order = c('TCRB training dataset', 'TCRB training dataset genes\n(not included during model training)', 'TCRB testing dataset', 'TCRA testing dataset', 'TCRG testing dataset', 'IGH testing dataset')
together$nice_loss_type = mapvalues(together$long_loss_type, from = loss_types, to=nice_loss_types)
together$nice_loss_type = factor(together$nice_loss_type, levels = loss_order)

# make trimming type names neat
together$trim_type = mapvalues(together$trim_type, from = unique(together$trim_type), to=c('J-gene trimming', 'V-gene trimming'))
together$trim_type = factor(together$trim_type, levels = c('V-gene trimming', 'J-gene trimming'))

# set plot path
ANNOTATION_TYPE <<- 'igor'
path = get_manuscript_path()

# set plot colors
require(RColorBrewer)
colors = brewer.pal(8, 'Set1')
colors = colors[!(colors %in% c("#FFFF33"))]
names(colors) = unique(nice_loss_types) 

# create plot to compare rel importance
plot4 = ggplot(together) +
    geom_hline(yintercept = 1, size = 3, color = 'gray60', linetype = 'dashed') +
    geom_vline(xintercept = 1, size = 3, color = 'gray60', linetype = 'dashed') +
    geom_point(aes(x = base_count_score, y = motif_score, color = nice_loss_type, shape = productivity), size = 10) +
    facet_grid(cols = vars(trim_type))+
    theme_cowplot(font_family = 'Arial') + 
    xlab('\nScale coefficient for base-count-beyond terms') +
    ylab('Scale coefficient for motif terms\n')+ 
    labs(color = 'Dataset', shape = 'Productivity') +
    background_grid(major = 'xy') + 
    scale_color_manual(values = colors, breaks = unique(nice_loss_types))+
    panel_border(color = 'gray60', size = 1.1) +
    theme(text = element_text(size = 38), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 38), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), strip.text = element_text(size = 38), panel.spacing = unit(3, "lines"), legend.key.height = unit(2, 'cm'), legend.key.width = unit(1.5, 'cm'))  

# get variable for training data set loss
training_loss = together[loss_type == 'igor' & productivity == 'nonproductive' & trim_type == 'V-gene trimming']$loss

# create plot to compare losses
plot5 = ggplot(together) +
    geom_hline(yintercept = training_loss, size = 3, color = 'gray60', linetype = 'dashed') +
    geom_vline(xintercept = training_loss, size = 3, color = 'gray60', linetype = 'dashed') +
    geom_abline(intercept = 0, size = 4, color = 'black')+
    geom_point(aes(x = loss, y = new_model_loss, color = nice_loss_type, shape = productivity), size = 10) +
    facet_grid(cols = vars(trim_type))+
    theme_cowplot(font_family = 'Arial') + 
    xlab('\nExpected per-sequence log loss (original model)') +
    ylab('Expected per-sequence log loss (scaled model)\n')+ 
    labs(color = 'Dataset', shape = 'Productivity') +
    background_grid(major = 'xy') + 
    scale_color_manual(values = colors, breaks = unique(nice_loss_types))+
    panel_border(color = 'gray60', size = 1.1) +
    theme(text = element_text(size = 38), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 38), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), strip.text = element_text(size = 38), panel.spacing = unit(3, "lines"), legend.key.height = unit(2, 'cm'), legend.key.width = unit(1.5, 'cm'))  

# isolate legend
legend = get_legend(plot4)
plot4f = plot4 + theme(legend.position = 'none')
plot5f = plot5 + theme(legend.position = 'none')

# align and combine plots
all = align_plots(plot4f, plot5f, align = 'vh', axis = 'lb')

first_grid = plot_grid(all[[1]], NULL, all[[2]], nrow = 3, labels = c('A', '', 'B'), label_size = 45, rel_heights = c(1, 0.03, 1)) 

tog = plot_grid(first_grid, NULL, legend, ncol = 3, rel_widths = c(1, 0.03, 0.4))

# save plot
file_name = paste0(path, '/feature_scale_together.pdf')
ggsave(file_name, plot = tog, width = 30, height = 27, units = 'in', dpi = 750, device = cairo_pdf)
