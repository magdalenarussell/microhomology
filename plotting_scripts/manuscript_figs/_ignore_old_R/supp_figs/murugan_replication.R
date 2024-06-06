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

args = commandArgs(trailingOnly=TRUE)

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
LEFT_NUC_MOTIF_COUNT <<- as.numeric(2)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(4)

UPPER_TRIM_BOUND <<- as.numeric(14) 

MODEL_TYPE <<- 'motif' 

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA

source(paste0(MOD_PROJECT_PATH,'scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'scripts/model_fitting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'plotting_scripts/plotting_functions.R'))

# get motif data
motif_data = aggregate_all_subject_data()
m_genes = c('TRBV2', 'TRBV7-2*01', 'TRBV10-3')

# Read in dist data
predicted_trims = get_predicted_distribution_data() 

# plot distribution for each gene
for (gene in unique(m_genes)){
    index = which(m_genes == gene)
    temp_plot = plot_predicted_trimming_dists(predicted_trims, gene, ylim = 0.75, color = '#01665e')
    temp_plot = temp_plot + 
        ggtitle(gene) +
        xlab('') +
        ylab('') +
        ylim(0, 0.8)
    assign(paste0('gene', index), temp_plot)
}

# Read in model coefficient data 
model_coef_matrix = get_model_coefficient_data() 
model_coef_matrix = model_coef_matrix[parameter %like% 'motif']

# set plotting positions
position_values = map_positions_to_values(unique(model_coef_matrix$parameter))
together = merge(model_coef_matrix, position_values, by.x = 'parameter', by.y = 'positions')

# convert to log_10 scale
together$log_10_pdel = together$coefficient/log(10)

# order variables
together$base = factor(together$base, levels = c('T', 'G', 'C', 'A'))
together$values = factor(together$values, levels = c(seq(4, 1), -1, -2))

# set bounds
limits = c(-0.516, 0.516)
motif_length = RIGHT_NUC_MOTIF_COUNT + LEFT_NUC_MOTIF_COUNT

# create heatmap
plot = ggplot(together, aes(x=values, y=base, fill=log_10_pdel)) +
    geom_tile() +
    theme_cowplot(font_family = 'Arial') + 
    xlab('Position') +
    ylab ('Base') +
    theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank()) +
    geom_vline(xintercept = RIGHT_NUC_MOTIF_COUNT + 0.5, size = 3.5, color = 'black') +
    guides(fill = guide_colourbar(barheight = 14)) +
    scale_fill_distiller(palette = 'PuOr', name = 'log10(probability of deletion)', limits = limits) +
    annotate("text", x = motif_length + 0.65, y = 0.25, label = "5\'", size = 8) +  
    annotate("text", x = 0.35, y = 0.25, label = "3\'", size = 8) +  
    coord_cartesian(ylim = c(1, 4), clip = "off")
 
heatmap = plot + 
    theme(text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 20)) 

# align and combine plots
all = align_plots(gene1, gene2, gene3, heatmap, align = 'v', axis = 'l')
first_grid = plot_grid(all[[1]], all[[2]], all[[3]], nrow = 1) +
    draw_label("Number of trimmed nucleotides", x=0.5, y=  0, vjust=0 , angle= 0, size = 30, fontfamily = 'Arial') +
    draw_label("Probability", x=  0, y=0.5, vjust= 0.5, angle=90, size = 30, fontfamily = 'Arial')

temp_grid = plot_grid(NULL, first_grid, ncol = 2, rel_widths = c(0.017, 1))
second_grid = plot_grid(temp_grid, NULL, all[[4]], nrow = 3, labels = c("A", "", "B"), label_size = 35, rel_heights = c(1, 0.05, 1.2), align = 'v')

# save plot
path = get_manuscript_path()
file_name = paste0(path, '/murugan_compare.pdf')
ggsave(file_name, plot = second_grid, width = 20, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
