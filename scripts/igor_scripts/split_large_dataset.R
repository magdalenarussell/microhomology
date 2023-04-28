library(data.table)
library(tidyverse)

source('config/config.R')
source('scripts/igor_scripts/igor_processing_functions.R')

args = commandArgs(trailingOnly=TRUE)

file = args[1]
output_directory = args[2]
size = as.numeric(args[3])

subject_ID = extract_subject_ID(file)
data = fread(file)
nonprod_subset = data[productive == FALSE]
total = nrow(nonprod_subset)

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

for (iter in seq(1, total, by = size)){
    end = min(iter + size-1, total)
    subset = nonprod_subset[iter:end, ]

    fwrite(subset, file.path(output_directory, paste0(subject_ID, '_', iter, '.tsv')), sep = '\t')
    print(paste0('finished iteration ', iter))
}

