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

MODEL_TYPE <<- 'motif_two-side-base-count-beyond'

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(args[11])

source(paste0(MOD_PROJECT_PATH, '/scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/model_fitting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

# get motif data
motif_data = aggregate_all_subject_data()
ordered_genes = motif_data[, .N, by = .(gene, p_gene)][order(-p_gene)]
top_genes = ordered_genes[1:6]$gene

# Read in dist data
predicted_trims = get_predicted_distribution_data() 

# plot a distribution for each gene
for (gene in unique(top_genes)){
    index = which(top_genes == gene)
    temp_plot = plot_predicted_trimming_dists(predicted_trims, gene, ylim = 0.95, color = '#01665e', motif_highlight_color = '#E79737', motif_highlight = c('CTT'))
    temp_plot = temp_plot + 
        ggtitle(gene) +
        xlab('') +
        ylab('') +
        ylim(0, 1)
    assign(paste0('gene', index), temp_plot)
}

# Read in model coefficient data 
pwm = get_model_coefficient_data() 

# plot heatmap of coefficients
heatmap = plot_model_coefficient_heatmap_single_group(pwm, with_values = FALSE, write_plot = FALSE, limits = c(-0.3, 0.3))
heatmap = heatmap + 
    theme(legend.position = 'none', text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 20)) 

heatmap2 = plot_base_count_coefficient_heatmap_single_group(pwm, with_values = FALSE, write_plot = FALSE, limits = c(-0.3, 0.3))
heatmap2 = heatmap2 + 
    theme(text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 20), legend.position = 'bottom', legend.direction = 'horizontal', legend.justification="center") +
    guides(fill = guide_colourbar(barwidth = 35, barheight = 2))

# isolate legend
legend = get_legend(heatmap2) 
heatmap2 = heatmap2 + theme(legend.position = 'none')

# align and combine all plots
all = align_plots(gene1, gene2, gene3, gene4, gene5, gene6, heatmap, heatmap2, legend, align = 'vh', axis = 'lbr')
first_grid = plot_grid(all[[1]], all[[2]], all[[3]], all[[4]], all[[5]], all[[6]], nrow = 2) +
    draw_label("Number of trimmed nucleotides", x=0.5, y=  0, vjust=0 , angle= 0, size = 30, fontfamily = 'Arial') +
    draw_label("Probability", x=  0, y=0.5, vjust= 0.5, angle=90, size = 30, fontfamily = 'Arial')
temp_grid = plot_grid(NULL, first_grid, ncol = 2, rel_widths = c(0.017, 1))
second_grid = plot_grid(all[[7]], all[[8]], nrow = 1, labels = c("B", "C"), label_size = 35, rel_widths = c(1, 1), align = 'h')
second_grid = plot_grid(second_grid, NULL, legend, nrow = 3, rel_heights = c(1, 0.02,0.2)) 
together = plot_grid(temp_grid, NULL, second_grid, nrow = 3, labels = c("A", "", ""), label_size = 35, rel_heights = c(0.7, 0.05, 0.6), align = 'v')

# save plot
path = get_manuscript_path()
file_name = paste0(path, '/motif_two-side-base-count.pdf')
ggsave(file_name, plot = together, width = 18, height = 13, units = 'in', dpi = 750, device = cairo_pdf)
