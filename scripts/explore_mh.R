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
LOWER_TRIM_BOUND <<- 2
UPPER_TRIM_BOUND <<- 14

blas_set_num_threads(NCPU)

path = paste0(PROJECT_PATH, '/plots/trim_by_join/',PRODUCTIVITY, '/', LOCUS, '/', TRIM_TYPE)
dir.create(path, recursive = TRUE)

source(paste0(PROJECT_PATH, '/scripts/data_processing_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/joining_gene_similarity_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/trimming_distribution_similarity_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/plotting_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/mh_functions.R'))

# get data
rep_data = read_all_data(directory = DATA_DIR)
rep_data = convert_adaptive_style_to_imgt(rep_data) 
rep_data_subset = filter_data(rep_data, filter_frequency = TRUE)
subset = rep_data_subset[v_trim >= LOWER_TRIM_BOUND & v_trim <= UPPER_TRIM_BOUND & j_trim >= LOWER_TRIM_BOUND & j_trim <= UPPER_TRIM_BOUND & n1_insertions == 0]
c_subset = subset[, .N, by = .(v_gene, j_gene, v_trim, j_trim, n1_insertions, productive)]

# get mh
c_subset = get_possible_mh(c_subset) 
c_subset = count_mh_bordering_trim(c_subset)
