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
JOINING_GENE <<- args[3]
DATA_DIR <<- args[4]
PRODUCTIVITY <<- 'nonproductive' 
NCPU <<- as.numeric(10)
GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')

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

# order genes by frequency
top_genes = unique(condensed_rep_data[order(-avg_paired_freq)][[GENE_NAME]])

# Create a trimming distribution plot for each gene
all_plots = c()
multinom_result = data.table()
model_fits = list()
model_coefs = data.table()
model_predictions = data.table()
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

    models = fit_multinom_models(uncondensed, joining_gene = 'j_gene') 

    result = get_pval(models$null, models$model)
    result$gene = gene_name

    coefs = convert_model_coefs_to_dt(models$model) 
    coefs$gene = gene_name
    model_coefs = rbind(model_coefs, coefs)

    preds = get_predicted_probs(models$model, uncondensed)
    model_predictions = rbind(model_predictions, preds)

    model_fits[[gene_name]] = models$model

    multinom_result = rbind(multinom_result, result)
    result = multinom_result[gene == gene_name]
    max_x = max(rep_data_subset[[TRIM_TYPE]]) - 5
    temp_plot = ggplot() +
        geom_line(data = empirical_data, aes(x = get(TRIM_TYPE), y = p_trim_given_pair, group = get(JOINING_GENE)), size = 0.6, alpha = 0.5) +
        geom_vline(xintercept = 0, color = 'black', size = 2) +
        geom_text(data = result, x = max_x, y = 0.15, aes(label = paste0('p = ', signif(p, 4))), size = 9) +
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

# write results
name = get_multinom_file_name(type = 'lrt_pvalue')
fwrite(multinom_result, name, sep = '\t')

model_name = get_multinom_file_name(type = 'model_coefs')
fwrite(model_coefs, model_name, sep = '\t')

prediction_name = get_multinom_file_name(type = 'model_preds')
fwrite(model_predictions, prediction_name, sep = '\t')

# combine all distributions
all = plot_grid(plotlist=mget(paste0("gene", index_list)), ncol = 5) 

# save plot
file_name = paste0(path, '/', LOCUS, '_', TRIM_TYPE, '_gallery_by_joining_', JOINING_GENE, '_signif.pdf')
ggsave(file_name, plot = all, width = 35, height = 30, units = 'in', dpi = 750, device = cairo_pdf, limitsize = FALSE)
