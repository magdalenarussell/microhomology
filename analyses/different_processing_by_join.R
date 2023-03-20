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
if (TRIM_TYPE %like% 'trim'){
    GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
} else {
    type = c('v_gene', 'j_gene')
    GENE_NAME <<- type[type != JOINING_GENE] 
}

blas_set_num_threads(NCPU)

path = paste0(PROJECT_PATH, '/plots/trim_by_join/',PRODUCTIVITY, '/', LOCUS, '/', TRIM_TYPE)
dir.create(path, recursive = TRUE)

source(paste0(PROJECT_PATH, '/scripts/data_processing_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/joining_gene_similarity_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/trimming_distribution_similarity_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/plotting_functions.R'))

# get data
rep_data = read_all_data(directory = DATA_DIR)
rep_data = convert_adaptive_style_to_imgt(rep_data) 
rep_data_subset = filter_data(rep_data, filter_frequency = TRUE)
preds = condense_data(rep_data_subset, filter = TRUE) 

# get pairwise trim dist differences
pairwise = pairwise_diffs(preds) 

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
 
file_name = paste0(path, '/pairwise_SAD_versus_hamming_', NT_COUNT, '_nt_', JOINING_GENE, '.pdf')
ggsave(file_name, plot = pairwise_dists_hamming, width = 35, height = 40, units = 'in', dpi = 750, device = cairo_pdf, limitsize = FALSE)

# plot pairwise trim dists by pairwise hamming
if (NT_COUNT > 5){
    dists_hamming_a = merge(pairwise, pairwise_hamming_aligned)

    pairwise_dists_hamming_a = plot_general_scatter(dists_hamming_a, yvar = 'sum_abs_diff', xvar = 'dist', ytitle = '\nPairwise sum of absolute differences in trimming distribution', xtitle = paste0('Pairwise hamming distance (first ', NT_COUNT, ' nt) using aligned sequences\n'), title = '', facet_var = GENE_NAME, facet_col = 5, add_trend = TRUE)
     
    file_name = paste0(path, '/pairwise_SAD_versus_hamming_aligned_', NT_COUNT, '_nt_', JOINING_GENE, '.pdf')
    ggsave(file_name, plot = pairwise_dists_hamming_a, width = 35, height = 40, units = 'in', dpi = 750, device = cairo_pdf, limitsize = FALSE)
}

# plot pairwise trim dists by pairwise alignment 
dists_align = merge(pairwise, pairwise_align)

pairwise_dists_align = plot_general_scatter(dists_align, yvar = 'sum_abs_diff', xvar = 'pairwise_score', ytitle = '\nPairwise sum of absolute differences in trimming distribution', xtitle = paste0('Pairwise global alignment score (first ', NT_COUNT, ' nt)\n'), title = '', facet_var = GENE_NAME, facet_col = 5, add_trend = TRUE)
 
file_name = paste0(path, '/pairwise_SAD_versus_alignment_', NT_COUNT, '_nt_', JOINING_GENE, '.pdf')
ggsave(file_name, plot = pairwise_dists_align, width = 35, height = 40, units = 'in', dpi = 750, device = cairo_pdf, limitsize = FALSE)
