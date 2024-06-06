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

ANNOTATION_TYPE<<- 'igor' 
TRIM_TYPE <<- 'v_trim' 
PRODUCTIVITY <<- 'nonproductive' 
MOTIF_TYPE <<- 'unbounded' 
NCPU <<- as.numeric(2)
GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
MODEL_GROUP <<- 'all_subjects' 
GENE_WEIGHT_TYPE <<- 'p_gene_marginal' 

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 0 
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 0

UPPER_TRIM_BOUND <<- as.numeric(14) 
LOWER_TRIM_BOUND <<- 2 

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA

MODEL_TYPE <<- 'distance'

types = c('v_gene_family_loss', 'full_v_gene_family_loss') 

# create a gene tree for each "most different" sequence protocol
for (type in types) {
    TYPE <<- type
    source(paste0(MOD_PROJECT_PATH, 'scripts/data_compilation_functions.R'))
    source(paste0(MOD_PROJECT_PATH,'scripts/model_fitting_functions.R'))
    source(paste0(MOD_PROJECT_PATH,'plotting_scripts/plotting_functions.R'))
    source(paste0(MOD_PROJECT_PATH,'scripts/model_evaluation_functions.R'))

    # assign cluster count and "most different" sequence protocol specific parameters
    if (TYPE == 'v_gene_family_loss'){
        cluster_count = 2
        combine_by_terminal = TRUE
        full_sequence = FALSE
        align = FALSE
    } else if (TYPE == 'full_v_gene_family_loss'){
        cluster_count = 4 
        combine_by_terminal = FALSE
        full_sequence = TRUE
        align = TRUE 
    }
    
    # gene the gene families and clusters
    v_families = get_gene_families(cluster_count, combine_by_terminal, full_sequence, align)
    clusters_grouped = v_families$cluster_data$clusters_grouped
    clusters_grouped[clusters_grouped != 1] = 2
    names(clusters_grouped) = v_families$cluster_data$gene
    
    # create tree
    require(RColorBrewer)
    colors = brewer.pal(cluster_count, 'Set2')

    plot(v_families$tree, type = 'unrooted', tip.color = colors[clusters_grouped], no.margin = TRUE) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    temp_plot = recordPlot()
    plot.new()
    assign(type, temp_plot)
}

# combine trees
together = plot_grid(v_gene_family_loss, full_v_gene_family_loss, nrow = 1, labels = c("A", "B"), label_size = 35, rel_heights = c(1, 1), align = 'h')

# save plot
path = get_manuscript_path()
file_name = paste0(path, '/gene_trees.pdf')
ggsave(file_name, plot = together, width = 20, height = 8, units = 'in', dpi = 750, device = cairo_pdf)


