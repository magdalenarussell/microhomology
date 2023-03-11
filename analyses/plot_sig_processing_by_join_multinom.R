source('config/config.R')

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(cowplot)
library(RhpcBLASctl)

args = commandArgs(trailingOnly=TRUE)

LOCUS <<- args[1]
stopifnot(LOCUS %in% c('TRB', 'TRA', 'TRA_igor'))
TRIM_TYPE <<- args[2] 
stopifnot(!(TRIM_TYPE %like% 'trim'))
JOINING_GENE <<- args[3]
DATA_DIR <<- args[4]
PRODUCTIVITY <<- 'nonproductive' 
NCPU <<- as.numeric(10)
if (TRIM_TYPE %like% 'trim'){
    GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
} else {
    type = c('v_gene', 'j_gene')
    GENE_NAME <<- type[type != JOINING_GENE] 
}

blas_set_num_threads(NCPU)

path = paste0(PROJECT_PATH, '/plots/trim_by_join')
dir.create(path, recursive = TRUE)

source(paste0(PROJECT_PATH, '/scripts/data_processing_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/joining_gene_similarity_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/plotting_functions.R'))

# Compile repertoire data for all subjects
rep_data = read_all_data(directory = DATA_DIR)
rep_data = convert_adaptive_style_to_imgt(rep_data) 
rep_data_subset = filter_data(rep_data, filter_frequency = TRUE)
condensed_rep_data = condense_data(rep_data_subset, filter = TRUE) 
true_condensed = condensed_rep_data[get(TRIM_TYPE) == 'TRUE']

# get results
name = get_multinom_file_name(type = 'lrt_pvalue')
multinom_result = fread(name)

top_genes = multinom_result[order(p)]$gene

# Create a trimming distribution plot for each gene
all_plots = c()
index_list = c()

for (gene_name in unique(top_genes)){
    print(paste0('starting ', gene_name))
    index = which(top_genes == gene_name)
    empirical_data = true_condensed[get(GENE_NAME) == gene_name]
    uncondensed = rep_data_subset[get(GENE_NAME) == gene_name]
    
    subset = unique(empirical_data[[JOINING_GENE]])

    if (length(subset) > 1){
        index_list = c(index_list, index)
    } else {
        next
    }

    result = multinom_result[gene == gene_name]

    xlabel = get_xaxis_multinom()

    max_x = median(empirical_data$avg_prob)

    temp_plot = ggplot() +
        geom_histogram(data = empirical_data, aes(x = avg_prob), alpha = 0.6) +
        geom_text(data = result, x = max_x, y = 0.15, aes(label = paste0('p = ', signif(p, 4))), size = 9) +
        ggtitle(gene_name) +
        xlab(xlabel) +
        ylab('Probability') +
        theme_cowplot(font_family = 'Arial') + 
    theme(text = element_text(size = 25), axis.text.x=element_text(size = 20), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5), legend.position = 'bottom', legend.direction = 'horizontal', legend.justification="center") + 
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5)
 
    legend = get_legend(temp_plot)
    temp_plot = temp_plot + 
        xlab('') +
        ylab('') +
        theme(legend.position = "none", text = element_text(size = 20))

    assign(paste0('gene', index), temp_plot)
}

# combine all distributions
all = plot_grid(plotlist=mget(paste0("gene", index_list)), ncol = 5) 

# save plot
file_name = paste0(path, '/', PRODUCTIVITY, '/', LOCUS, '_', TRIM_TYPE, '_gallery_by_joining_', JOINING_GENE, '_signif.pdf')
ggsave(file_name, plot = all, width = 35, height = 30, units = 'in', dpi = 750, device = cairo_pdf, limitsize = FALSE)
