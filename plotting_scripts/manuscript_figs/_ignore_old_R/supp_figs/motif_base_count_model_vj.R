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
TRIM_TYPE <<- 'j_trim' 
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

# get motif data
j_motif_data = aggregate_all_subject_data()
j_ordered_genes = j_motif_data[, .N, by = .(gene, p_gene)][order(-p_gene)]
j_top_genes = j_ordered_genes[1:6]$gene

# Read in dist data
j_predicted_trims = get_predicted_distribution_data() 

# plot j-gene trimming distributions
for (gene in unique(j_top_genes)){
    index = which(j_top_genes == gene)
    j_temp_plot = plot_predicted_trimming_dists(j_predicted_trims, gene, ylim = 0.75, color = '#01665e')
    j_temp_plot = j_temp_plot + 
        ggtitle(gene) +
        xlab('') +
        ylab('') +
        ylim(0, 0.8)
    assign(paste0('jgene', index), j_temp_plot)
}

# Read in model coefficient data 
j_pwm = get_model_coefficient_data() 

# plot a heatmap for model coefficients
j_heatmap = plot_model_coefficient_heatmap_single_group(j_pwm, with_values = FALSE, write_plot = FALSE, limits = c(-0.24, 0.24))
j_heatmap = j_heatmap + 
    theme(legend.position = 'none', text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 20)) 

j_heatmap2 = plot_base_count_coefficient_heatmap_single_group(j_pwm, with_values = FALSE, write_plot = FALSE, limits = c(-0.24, 0.24))
j_heatmap2 = j_heatmap2 + 
    theme(text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 20)) 

# align and combine all plots
all = align_plots(jgene1, jgene2, jgene3, jgene4, jgene5, jgene6, j_heatmap, j_heatmap2, align = 'vh', axis = 'lb')

first_grid = plot_grid(all[[1]], all[[2]], all[[3]], all[[4]], all[[5]], all[[6]], nrow = 2) +
    draw_label("Number of trimmed nucleotides", x=0.5, y=  0, vjust=0 , angle= 0, size = 30, fontfamily = 'Arial') +
    draw_label("Probability", x=  0, y=0.5, vjust= 0.5, angle=90, size = 30, fontfamily = 'Arial')

temp_grid = plot_grid(NULL, first_grid, ncol = 2, rel_widths = c(0.017, 1))
second_grid = plot_grid(all[[7]], all[[8]], nrow = 1, labels = c("B", "C"), label_size = 35, rel_widths = c(0.55, 1), align = 'h')

together = plot_grid(temp_grid, NULL, second_grid, nrow = 3, labels = c("A", "", ""), label_size = 35, rel_heights = c(1, 0.05, 0.6), align = 'v')

# save plot
TRIM_TYPE <<- 'v_trim'
path = get_manuscript_path()
file_name = paste0(path, '/motif_two-side-base-count_j.pdf')
ggsave(file_name, plot = together, width = 20, height = 15, units = 'in', dpi = 750, device = cairo_pdf)
