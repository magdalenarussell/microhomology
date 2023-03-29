library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)

args = commandArgs(trailingOnly=TRUE)

output_directory = args[1]
NCPU <<- args[2]

source('config/config.R')
source('scripts/igor_scripts/igor_processing_functions.R')

output = fread(file.path(output_directory, 'foo_output/best_scenarios_counts.csv'))

sample = sample_sequences(output)

fwrite(sample, file.path(output_directory, 'foo_output/sampled_scenarios.csv'), sep = ';')



