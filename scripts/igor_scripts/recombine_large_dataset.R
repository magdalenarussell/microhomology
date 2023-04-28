library(data.table)
library(tidyverse)
library(vroom)

source('config/config.R')
source('scripts/igor_scripts/igor_processing_functions.R')

args = commandArgs(trailingOnly=TRUE)

input_directory = args[1]
output_directory = args[2]

files = list.files(input_directory, full.names = TRUE)
subject_ID = extract_subject_ID(files[1])
data = vroom(files)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

fwrite(data, file.path(output_directory, paste0(subject_ID, '.tsv')), sep ='\t')
