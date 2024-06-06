library(data.table)
# library(tidyverse)

source('config/config.R')
source('scripts/igor_processing_functions.R')

args = commandArgs(trailingOnly=TRUE)

file = args[1]
output_directory = args[2]
adaptive = as.logical(args[3])
LOCUS <<- args[4]
NCPU <<- as.numeric(args[5])
stopifnot(LOCUS %in% c('alpha', 'beta'))

subject_ID = extract_subject_ID(file)
cdr3_data = get_raw_cdr3_seqs(file, adaptive)
cdr3s = cdr3_data$cdr3s
v_index = cdr3_data$v_index
counts = cdr3_data$counts

output_location = file.path(output_directory, subject_ID)
dir.create(output_location, recursive = TRUE, showWarnings = FALSE)

fwrite(as.data.table(cdr3s), file.path(output_location, 'cdr3_seqs.txt'), sep = '\t', col.names = FALSE)
write.csv2(as.data.frame(as.numeric(v_index)), file.path(output_location, 'v_anchor.csv'), col.names = FALSE, quote = FALSE, row.names = FALSE)
fwrite(as.data.table(counts), file.path(output_location, 'cdr3_counts.txt'), sep = '\t', col.names = FALSE)

cat(output_location)
