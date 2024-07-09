source('config/config.R')
source('config/file_paths.R')

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(ggplot2)
library(cowplot)
library(Biostrings)
library(mclogit)
library(matrixcalc)
library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)

args = commandArgs(trailingOnly=TRUE)

ANNOTATION_TYPE <<- 'igor_alpha'
PARAM_GROUP <<- args[1]
PARAM_GROUP <<- 'nonproductive_v-j_trim_ligation-mh'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))
NCPU<<-2
# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 1
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 2
MODEL_TYPE <<- 'motif_two-side-base-count-beyond_average-mh_ligation-mh'
L2 <<- 'True'

source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

predictions_file = get_validation_predictions_file_path(L2, validation_annotation = 'igor_prob_experiment', sample_annotation=SAMPLE_ANNOT)

pred = fread(predictions_file)[!(new_type == '0')]

MODEL_TYPE <<- 'motif_two-side-base-count-beyond'
source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

predictions_file_noMH = get_validation_predictions_file_path(L2, validation_annotation = 'igor_prob_experiment', sample_annotation=SAMPLE_ANNOT)

pred_noMH = fread(predictions_file_noMH)[!(new_type == '0')]
setnames(pred_noMH, 'predicted_prob', 'predicted_prob_noMH')

# merge
tog = merge(pred, pred_noMH)

tog[, log_diff := log(predicted_prob) - log(predicted_prob_noMH)]
tog[, diff := predicted_prob - predicted_prob_noMH]

# subset by gene usage
# get gene probabilities from IGoR (filter out low prob genes)
j_usage = fread(paste0(MOD_OUTPUT_PATH, '/igor_alpha/nonproductive_v-j_trim_ligation-mh/unbounded_motif_trims_bounded_-2_14/1_2_motif_two-side-base-count-beyond_average-mh_ligation-mh/igor_prob_experiment/igor_alpha_jchoice_params.tsv'))
v_usage = fread(paste0(MOD_OUTPUT_PATH, '/igor_alpha/nonproductive_v-j_trim_ligation-mh/unbounded_motif_trims_bounded_-2_14/1_2_motif_two-side-base-count-beyond_average-mh_ligation-mh/igor_prob_experiment/igor_alpha_vchoice_params.tsv'))
j_usage$indicator = 1
v_usage$indicator = 1

joint_usage = merge(v_usage, j_usage, by = 'indicator', allow.cartesian = TRUE)
joint_usage[, joint_gene_usage := v_gene_prob * j_gene_prob]

## merge and restrict to genes that have high usage in IGoR
tog_merge = merge(tog, joint_usage, by = c('v_gene', 'j_gene'))
tog_subset = tog_merge[joint_gene_usage > 1e-4]

plot_trend = ggplot(tog_subset) +
    geom_hex(aes(x = ligation_mh, y = diff)) +
    geom_abline(intercept = 0, slope = 0, color = 'black', linewidth = 1) +
    geom_smooth(aes(x = ligation_mh, y = diff), method = 'lm', se = TRUE, linewidth = 1) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('MH within trimming/ligation annotation') +
    ylab('MH model probability - noMH model probability\n(MH annotation space)') +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    scale_fill_gradient(low = '#efedf5', high = '#3f007d', trans = 'log10', name = expression(log[10]("count"))) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_igor_comparison/TRA_mh_comparison_trend_ligationMH.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = plot_trend, width = 10, height = 7, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)
