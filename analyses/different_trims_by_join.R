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
NT_COUNT <<- args[5]
if (NT_COUNT != 'all'){
    NT_COUNT <<- as.numeric(NT_COUNT)
}

PRODUCTIVITY <<- args[6] 

NCPU <<- as.numeric(10)
GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')

blas_set_num_threads(NCPU)

path = paste0(PROJECT_PATH, '/plots/trim_by_join')
dir.create(path, recursive = TRUE)

source(paste0(PROJECT_PATH, '/scripts/data_processing_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/joining_gene_similarity_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/trimming_distribution_similarity_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/plotting_functions.R'))

file = prediction_name = get_multinom_file_name(type = 'model_preds') 
preds = fread(file)

# get pairwise trim dist differences
pairwise = pairwise_diffs(preds) 

# cluster trim dists
cluster = cluster_diffs(preds, cluster_count = 2)
pairwise_cluster = transform_cluster_data_to_pairwise(cluster)

# compare to sequence-based features of the joining gene level
## get pairwise hamming
if (NT_COUNT > 5){
    pairwise_hamming_aligned = get_pairwise_hamming(joining_genes = unique(preds[[JOINING_GENE]]), nt_count = NT_COUNT, alignment = 'align')
}
pairwise_hamming = get_pairwise_hamming(joining_genes = unique(preds[[JOINING_GENE]]), nt_count = NT_COUNT, alignment = 'no_align')
pairwise_align = get_pairwise_msa(joining_genes = unique(preds[[JOINING_GENE]]), nt_count = NT_COUNT)


# plot pairwise trim dists by pairwise hamming
dists_hamming = merge(pairwise, pairwise_hamming)


pairwise_dists_hamming = plot_general_scatter(dists_hamming, yvar = 'sum_abs_diff', xvar = 'dist', ytitle = '\nPairwise sum of absolute differences in trimming distribution', xtitle = paste0('Pairwise hamming distance (first ', NT_COUNT, ' nt)\n'), title = '', facet_var = GENE_NAME, facet_col = 5, add_trend = TRUE)
 
file_name = paste0(path, '/', PRODUCTIVITY, '/', LOCUS, '_', TRIM_TYPE, '_pairwise_SAD_versus_hamming_', NT_COUNT, '_nt_', JOINING_GENE, '.pdf')

ggsave(file_name, plot = pairwise_dists_hamming, width = 35, height = 40, units = 'in', dpi = 750, device = cairo_pdf, limitsize = FALSE)

# plot pairwise trim dists by pairwise hamming
if (NT_COUNT > 5){
    dists_hamming_a = merge(pairwise, pairwise_hamming_aligned)

    pairwise_dists_hamming_a = plot_general_scatter(dists_hamming_a, yvar = 'sum_abs_diff', xvar = 'dist', ytitle = '\nPairwise sum of absolute differences in trimming distribution', xtitle = paste0('Pairwise hamming distance (first ', NT_COUNT, ' nt) using aligned sequences\n'), title = '', facet_var = GENE_NAME, facet_col = 5, add_trend = TRUE)
     
    file_name = paste0(path, '/', PRODUCTIVITY, '/', LOCUS, '_', TRIM_TYPE, '_pairwise_SAD_versus_hamming_aligned_', NT_COUNT, '_nt_', JOINING_GENE, '.pdf')

    ggsave(file_name, plot = pairwise_dists_hamming_a, width = 35, height = 40, units = 'in', dpi = 750, device = cairo_pdf, limitsize = FALSE)
}

# plot pairwise trim dists by pairwise alignment 
dists_align = merge(pairwise, pairwise_align)


pairwise_dists_align = plot_general_scatter(dists_align, yvar = 'sum_abs_diff', xvar = 'pairwise_score', ytitle = '\nPairwise sum of absolute differences in trimming distribution', xtitle = paste0('Pairwise global alignment score (first ', NT_COUNT, ' nt)\n'), title = '', facet_var = GENE_NAME, facet_col = 5, add_trend = TRUE)
 
file_name = paste0(path, '/', PRODUCTIVITY, '/', LOCUS, '_', TRIM_TYPE, '_pairwise_SAD_versus_alignment_', NT_COUNT, '_nt_', JOINING_GENE, '.pdf')

ggsave(file_name, plot = pairwise_dists_align, width = 35, height = 40, units = 'in', dpi = 750, device = cairo_pdf, limitsize = FALSE)


# plot clustered trim dists by pairwise hamming
cluster_hamming = merge(pairwise_cluster, pairwise_hamming)


cluster_dists_hamming = plot_general_boxplot(cluster_hamming[!(get(paste0(JOINING_GENE, '.x')) == get(paste0(JOINING_GENE, '.y')))], xvar = 'cluster', yvar = 'dist', xtitle = '\nTrimming distribution cluster (from K-means)', ytitle = paste0('Pairwise hamming distance (first ', NT_COUNT, ' nt)\n'), title = '', facet_var = GENE_NAME, facet_col = 5)
 
file_name = paste0(path, '/', PRODUCTIVITY, '/', LOCUS, '_', TRIM_TYPE, '_cluster_versus_hamming_', NT_COUNT, '_nt_', JOINING_GENE, '.pdf')

ggsave(file_name, plot = cluster_dists_hamming, width = 35, height = 40, units = 'in', dpi = 750, device = cairo_pdf, limitsize = FALSE)


# plot clustered trim dists by pairwise align 
cluster_align = merge(pairwise_cluster, pairwise_align)


cluster_dists_align = plot_general_boxplot(cluster_align[!(get(paste0(JOINING_GENE, '.x')) == get(paste0(JOINING_GENE, '.y')))], xvar = 'cluster', yvar = 'pairwise_score', xtitle = '\nTrimming distribution cluster (from K-means)', ytitle = paste0('Pairwise global alignment score (first ', NT_COUNT, ' nt)\n'), title = '', facet_var = GENE_NAME, facet_col = 5)
 
file_name = paste0(path, '/', PRODUCTIVITY, '/', LOCUS, '_', TRIM_TYPE, '_cluster_versus_align_', NT_COUNT, '_nt_', JOINING_GENE, '.pdf')

ggsave(file_name, plot = cluster_dists_align, width = 35, height = 40, units = 'in', dpi = 750, device = cairo_pdf, limitsize = FALSE)


