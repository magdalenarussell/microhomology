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

source(paste0(MOD_PROJECT_PATH, 'scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'scripts/model_fitting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'plotting_scripts/plotting_functions.R'))

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
    temp_plot = plot_predicted_trimming_dists(predicted_trims, gene, ylim = 0.75, color = '#01665e', motif_highlight_color = '#E79737', motif_highlight = c('CTT', 'CGT'), seq_text = 5)
    temp_plot = temp_plot + 
        ggtitle(gene) +
        xlab('') +
        ylab('') +
        ylim(0, 0.8)+
        theme(text = element_text(size = 20), plot.margin = margin(l = -0.1, b = -0.3, unit = 'cm'))
    assign(paste0('gene', index), temp_plot)
}

# combine all distributions
all = plot_grid(plotlist=mget(paste0("gene", seq(1, length(top_genes)))), ncol = 5)

# save plot
path = get_manuscript_path()
file_name = paste0(path, '/trb_gallery.pdf')
ggsave(file_name, plot = all, width = 30, height = 38, units = 'in', dpi = 750, device = cairo_pdf, limitsize = FALSE)


