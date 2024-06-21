library(data.table)
# library(tidyverse)
library(parallel)

source('config/config.R')
source('igor_annotation_scripts/igor_processing_functions.R')

args = commandArgs(trailingOnly=TRUE)

file = args[1]
output_directory = args[2]
adaptive = as.logical(args[3])
sample_size = args[4]
LOCUS <<- args[5]
NCPU <<- as.numeric(args[6])
stopifnot(LOCUS %in% c('alpha', 'beta'))

subject_ID = extract_subject_ID(file)
cdr3_data = get_raw_cdr3_seqs_first_subsample(file, sample_size, adaptive)
cdr3s = cdr3_data$cdr3s
v_index = cdr3_data$v_index
counts = cdr3_data$counts

output_location = file.path(output_directory)
dir.create(output_location, recursive = TRUE, showWarnings = FALSE)

fwrite(as.data.table(cdr3s), file.path(output_location, paste0(subject_ID, '_cdr3_seqs.txt')), sep = '\t', col.names = FALSE)

cat(output_location)
