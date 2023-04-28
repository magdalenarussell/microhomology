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

# get data
rep_data = read_all_data(directory = DATA_DIR)
rep_data_subset = process_data(rep_data)
subset = rep_data_subset[v_trim >= LOWER_TRIM_BOUND & v_trim <= UPPER_TRIM_BOUND & j_trim >= LOWER_TRIM_BOUND & j_trim <= UPPER_TRIM_BOUND & vj_insert == 0]
c_subset = subset[, .N, by = .(v_gene, j_gene, v_trim, j_trim, vj_insert, productive)]
c_subset[, pair := paste0(v_gene, '.', j_gene)]

germline = data.table()
for (pair in unique(c_subset$pair)){
    temp = data.table(v_trim = rep(seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND), each = UPPER_TRIM_BOUND-LOWER_TRIM_BOUND + 1), j_trim = rep(seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND), UPPER_TRIM_BOUND-LOWER_TRIM_BOUND + 1))
    temp2 = data.table()
    temp$pair = pair
    temp$v_gene = str_split(pair, '\\.')[[1]][1]
    temp$j_gene = str_split(pair, '\\.')[[1]][2]
    germline = rbind(germline, temp)
}

germline$productive = PRODUCTIVITY
germline$vj_insert = 0

tog = merge(germline, c_subset, fill = TRUE, all.x = TRUE)
tog[is.na(N), N := 0]

# get mh
tog = get_possible_mh(tog) 
tog = count_mh_bordering_trim(tog)

path = file.path('analyses', '_ignore')
dir.create(path)
fwrite(tog, paste0(path, '/', LOCUS, '_', TRIM_TYPE, '_', JOINING_GENE, '_bordering_mh_', DATA_TYPE, '_annotations.tsv'), sep = '\t') 
