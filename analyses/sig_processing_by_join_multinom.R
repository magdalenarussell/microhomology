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
stopifnot(LOCUS %in% c('TRA', 'TRA_igor'))
TRIM_TYPE <<- args[2] 
JOINING_GENE <<- args[3]
DATA_DIR <<- args[4]
PRODUCTIVITY <<- args[5] 
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

# order genes by frequency
top_genes = unique(condensed_rep_data[order(-avg_paired_freq)][[GENE_NAME]])

# Create a trimming distribution plot for each gene
multinom_result = data.table()
model_fits = list()
model_coefs = data.table()
model_predictions = data.table()

for (gene_name in unique(top_genes)){
    print(paste0('starting ', gene_name))
    empirical_data = condensed_rep_data[get(GENE_NAME) == gene_name]
    uncondensed = rep_data_subset[get(GENE_NAME) == gene_name]
    
    subset = unique(empirical_data[[JOINING_GENE]])

    if (!(length(subset) > 1)){
        next
    }

    models = fit_multinom_models(uncondensed, joining_gene = JOINING_GENE) 

    result = get_pval(models$null, models$model)
    result$gene = gene_name
    multinom_result = rbind(multinom_result, result)

    coefs = convert_model_coefs_to_dt(models$model) 
    coefs$gene = gene_name
    model_coefs = rbind(model_coefs, coefs)

    preds = get_predicted_probs(models$model, uncondensed)
    model_predictions = rbind(model_predictions, preds)

    model_fits[[gene_name]] = models$model
}

# write results
name = get_multinom_file_name(type = 'lrt_pvalue')
fwrite(multinom_result, name, sep = '\t')

model_name = get_multinom_file_name(type = 'model_coefs')
fwrite(model_coefs, model_name, sep = '\t')

prediction_name = get_multinom_file_name(type = 'model_preds')
fwrite(model_predictions, prediction_name, sep = '\t')
