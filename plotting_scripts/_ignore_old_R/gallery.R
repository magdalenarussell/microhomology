source('mechanistic-trimming/config/config.R')

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

MODEL_TYPE <<- args[11]

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(args[12])

source(paste0(MOD_PROJECT_PATH, '/scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/model_fitting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

# get motif data
motif_data = aggregate_all_subject_data()
ordered_genes = motif_data[p_gene > 1e-04, .N, by = .(gene, p_gene)][order(-p_gene)]
top_genes = ordered_genes$gene

# Read in dist data
predicted_trims = get_predicted_distribution_data() 

# Create a trimming distribution plot for each gene
all_plots = c()
for (gene in unique(top_genes)){
    index = which(top_genes == gene)
    temp_plot = plot_predicted_trimming_dists(predicted_trims, gene, ylim = 0.85, color = '#01665e', motif_highlight_color = '#E79737', motif_highlight = c('CTT', 'CGT'), seq_text = 5)
    temp_plot = temp_plot + 
        ggtitle(gene) +
        xlab('') +
        ylab('') +
        ylim(0, 0.9)+
        theme(text = element_text(size = 20), plot.margin = margin(l = -0.1, b = -0.3, unit = 'cm'))
    assign(paste0('gene', index), temp_plot)
}

# combine all distributions
all = plot_grid(plotlist=mget(paste0("gene", seq(1, length(top_genes)))), ncol = 5)

# save plot
path = get_manuscript_path()
file_name = paste0(path, '/tra_gallery.pdf')
ggsave(file_name, plot = all, width = 30, height = 38, units = 'in', dpi = 750, device = cairo_pdf, limitsize = FALSE)


