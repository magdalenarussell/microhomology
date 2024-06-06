source('mechanistic-trimming/config/config.R')

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(Biostrings)
library(ggplot2)
library(cowplot)
library(mclogit)
library(matrixcalc)
library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)

ANNOTATION_TYPE <<- 'igor_sim_alpha' 

PARAM_GROUP <<- 'nonproductive_v-j_trim_ligation-mh'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))

NCPU <<- as.numeric(2)

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(1)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(2)

MODEL_TYPE <<- 'motif_two-side-base-count-beyond_ligation-mh'

L2 <<- 'True'

source(paste0(MOD_PROJECT_PATH, '/scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

filename = processed_data_path()
processed = fread(filename)

ANNOTATION_TYPE <<- 'igor_sim_alpha_independent_vj_choice'

source(paste0(MOD_PROJECT_PATH, '/scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

filename = processed_data_path()
shuffled_genes = fread(filename)

ANNOTATION_TYPE <<- 'igor_sim_alpha_shuffled_trimming'

source(paste0(MOD_PROJECT_PATH, '/scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

filename = processed_data_path()
shuffled_trims = fread(filename)
shuffled_trims = shuffled_trims[, -23]

ANNOTATION_TYPE <<- 'igor_sim_random_sequence'

source(paste0(MOD_PROJECT_PATH, '/scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

filename = processed_data_path()
random_seq = fread(filename)
random_seq = random_seq[, -23]

ANNOTATION_TYPE <<- 'no_mh_model_sim_alpha'

source(paste0(MOD_PROJECT_PATH, '/scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

filename = processed_data_path()
no_mh_sim = fread(filename)
no_mh_sim = no_mh_sim[, -23]

# plot original with shuffled genes--V-genes
shuffled_genes_cond_v = shuffled_genes[, sum(count), by = .(v_gene, v_trim, ligation_mh)]
processed_cond_v = processed[, sum(count), by = .(v_gene, v_trim, ligation_mh)]
shuffled_genes_cond_v[, scenario_freq_shuffled_genes := V1/sum(V1), by = .(v_gene)]
processed_cond_v[, scenario_freq := V1/sum(V1), by = .(v_gene)]

tog_cond_v = merge(processed_cond_v, shuffled_genes_cond_v, by = c('v_gene', 'v_trim', 'ligation_mh'), all = TRUE)
tog_cond_v[is.na(tog_cond_v)] <- 0

plot = ggplot(tog_cond_v) +
    geom_point(aes(x = scenario_freq, y = scenario_freq_shuffled_genes, color = as.factor(ligation_mh)), size = 3, alpha = 0.5) +
    theme_classic()

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/igor_experiments/explore_signal/trimming_dists_shuffled_genes_v_trim.pdf')
ggsave(file_name, plot = plot, width = 10, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

# plot original with shuffled genes--J-genes
shuffled_genes_cond_j = shuffled_genes[, sum(count), by = .(j_gene, j_trim, ligation_mh)]
processed_cond_j = processed[, sum(count), by = .(j_gene, j_trim, ligation_mh)]
shuffled_genes_cond_j[, scenario_freq_shuffled_genes := V1/sum(V1), by = .(j_gene)]
processed_cond_j[, scenario_freq := V1/sum(V1), by = .(j_gene)]

tog_cond_j = merge(processed_cond_j, shuffled_genes_cond_j, by = c('j_gene', 'j_trim', 'ligation_mh'), all = TRUE)
tog_cond_j[is.na(tog_cond_j)] <- 0

plot = ggplot(tog_cond_j) +
    geom_point(aes(x = scenario_freq, y = scenario_freq_shuffled_genes, color = as.factor(ligation_mh)), size = 3, alpha = 0.5) +
    theme_classic()

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/igor_experiments/explore_signal/trimming_dists_shuffled_genes_j_trim.pdf')
ggsave(file_name, plot = plot, width = 10, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

# plot original with shuffled trims
shuffled_trims[, scenario_freq_shuffled_trims := count/sum(count), by = .(j_gene, v_gene)]
processed[, scenario_freq := count/sum(count), by = .(v_gene, j_gene)]

tog_cond = merge(processed, shuffled_trims, by = c('v_gene', 'j_gene', 'v_trim', 'j_trim', 'ligation_mh'), all = TRUE)
tog_cond[is.na(tog_cond)] <- 0

plot = ggplot(tog_cond) +
    geom_point(aes(x = scenario_freq, y = scenario_freq_shuffled_trims, color = as.factor(ligation_mh)), size = 3, alpha = 0.5) +
    theme_classic()

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/igor_experiments/explore_signal/trimming_dists_shuffled_trims.pdf')
ggsave(file_name, plot = plot, width = 10, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

plot = ggplot(tog_cond[ligation_mh > 0]) +
    geom_point(aes(x = scenario_freq, y = scenario_freq_shuffled_trims, color = as.factor(ligation_mh)), size = 3, alpha = 0.5) +
    theme_classic()

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/igor_experiments/explore_signal/trimming_dists_shuffled_trims_nonzero_mh.pdf')
ggsave(file_name, plot = plot, width = 10, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

# plot original with random_seq
random_seq[, scenario_freq_random_seq := count/sum(count), by = .(j_gene, v_gene)]
processed[, scenario_freq := count/sum(count), by = .(v_gene, j_gene)]

tog_cond = merge(processed, random_seq, by = c('v_gene', 'j_gene', 'v_trim', 'j_trim', 'ligation_mh'), all = TRUE)
tog_cond[is.na(tog_cond)] <- 0

plot = ggplot(tog_cond) +
    geom_point(aes(x = scenario_freq, y = scenario_freq_random_seq, color = as.factor(ligation_mh)), size = 3, alpha = 0.5) +
    theme_classic()

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/igor_experiments/explore_signal/trimming_dists_random_seq.pdf')
ggsave(file_name, plot = plot, width = 10, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

plot = ggplot(tog_cond[ligation_mh > 0]) +
    geom_point(aes(x = scenario_freq, y = scenario_freq_random_seq, color = as.factor(ligation_mh)), size = 3, alpha = 0.5) +
    theme_classic()

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/igor_experiments/explore_signal/trimming_dists_random_seq_nonzero_mh.pdf')
ggsave(file_name, plot = plot, width = 10, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

# plot trimming dists and ligationMH
no_mh_sim[, scenario_freq_no_mh_sim := count/sum(count), by = .(v_gene, j_gene)]

tog = merge(processed, no_mh_sim, by = c('v_gene', 'j_gene', 'v_trim', 'j_trim', 'ligation_mh'))
