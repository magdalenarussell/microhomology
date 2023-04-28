source('config/config.R')

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
library(cowplot)
library(RhpcBLASctl)

args = commandArgs(trailingOnly=TRUE)

DATA_TYPE <<- args[1]
LOCUS <<- args[2]
stopifnot(LOCUS %in% c('TRB', 'TRA', 'TRA_igor'))
TRIM_TYPE <<- args[3] 
JOINING_GENE <<- args[4]
DATA_DIR <<- args[5]
NT_COUNT <<- args[6]
if (NT_COUNT != 'all'){
    NT_COUNT <<- as.numeric(NT_COUNT)
}
PRODUCTIVITY <<- args[7] 
NCPU <<- as.numeric(args[8]) 
if (TRIM_TYPE %like% 'trim'){
    GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
} else {
    type = c('v_gene', 'j_gene')
    GENE_NAME <<- type[type != JOINING_GENE] 
}
LOWER_TRIM_BOUND <<- as.numeric(args[9]) 
UPPER_TRIM_BOUND <<- as.numeric(args[10]) 

blas_set_num_threads(NCPU)
setDTthreads(NCPU)

path = paste0(PROJECT_PATH, '/plots/trim_by_join/',PRODUCTIVITY, '/', LOCUS, '/', TRIM_TYPE)
dir.create(path, recursive = TRUE)

source(paste0(PROJECT_PATH, '/scripts/processing_functions/data_processing_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/analysis_functions/joining_gene_similarity_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/analysis_functions/trimming_distribution_similarity_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/plotting_functions/plotting_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/processing_functions/mh_functions.R'))

# Compile repertoire data for all subjects
rep_data = read_all_data(directory = DATA_DIR)
rep_data_subset = process_data(rep_data)
condensed_rep_data = condense_data(rep_data_subset, filter = TRUE) 

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
    empirical_data = condensed_rep_data[get(GENE_NAME) == gene_name]
    uncondensed = rep_data_subset[get(GENE_NAME) == gene_name]
    
    subset = unique(empirical_data[[JOINING_GENE]])

    if (length(subset) > 1){
        index_list = c(index_list, index)
    } else {
        next
    }

    result = multinom_result[gene == gene_name]

    temp_plot = ggplot() +
        geom_line(data = empirical_data, aes(x = get(TRIM_TYPE), y = avg_prob, group = get(JOINING_GENE)), size = 0.6, alpha = 0.5) +
        geom_vline(xintercept = 0, color = 'black', size = 2) +
        geom_text(data = result, x = Inf, y = Inf, aes(label = paste0('p = ', signif(p, 4))), size = 9, vjust = 1.3, hjust = 1.1) +
        ggtitle(gene_name) +
        xlab('Number of trimmed nucleotides') +
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
height = round(length(index_list)/5)*6.5

# save plot
file_name = paste0(path, '/', DATA_TYPE, '_gallery_by_joining_', JOINING_GENE, '_signif.pdf')
ggsave(file_name, plot = all, width = 35, height = height, units = 'in', dpi = 750, device = cairo_pdf, limitsize = FALSE)
