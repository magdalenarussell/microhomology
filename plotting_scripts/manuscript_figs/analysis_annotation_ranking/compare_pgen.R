source('config/config.R')
source('config/file_paths.R')

library(ggpubr)
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

ANNOTATION_TYPE <<- 'igor_alpha'
PARAM_GROUP <<- 'nonproductive_v-j_trim_ligation-mh'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))
NCPU <<- as.numeric(2)
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

noMH_predictions_file = get_validation_predictions_file_path(L2, validation_annotation = 'igor_prob_experiment', sample_annotation=SAMPLE_ANNOT)

no_pred = fread(noMH_predictions_file)[!(new_type == '0')]

# get average MH
pred[, avg_mh := mean(ligation_mh), by = .(index)]

# sum by sequence
pred_pgen = pred[, sum(predicted_prob), by = c('index', 'v_gene', 'j_gene', 'avg_mh')]
setnames(pred_pgen, 'V1', 'mh_pgen')

no_pred_pgen = no_pred[, sum(predicted_prob), by = c('index', 'v_gene', 'j_gene')]
setnames(no_pred_pgen, 'V1', 'no_mh_pgen')

# join everything together .x vars are mh model and .y vars are noMH model
tog = merge(pred_pgen, no_pred_pgen, by = c('index', 'v_gene', 'j_gene'))


compare = ggplot(tog, aes(y = log10(mh_pgen), x = log10(no_mh_pgen)))+
    geom_hex()+
    theme_cowplot(font_family = 'Arial') + 
    xlab('log10 sequence generation probability from no-MH model (per gene-pair)') +
    ylab('log10 sequence generation probability from MH model (per gene-pair)') +
    geom_abline(intercept = 0, slope = 1, color = 'black', linewidth = 1) +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 0.5) +
    scale_fill_gradient(low = '#e8f5f1', high = '#1b9e77', name = 'count')+
    theme(text = element_text(size = 16), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 12), plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), legend.position = 'bottom', legend.justification="center", legend.key.width = unit(1, "in"), legend.key.height = unit(1, "cm"), panel.grid.major = element_line(size = 0.25)) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_annotation_ranking/pgen_comparison.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = compare, width = 8.5, height = 8.5, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)

regress = lm(data = tog, log10(mh_pgen) ~ log10(no_mh_pgen))
print(regress)

# Paired t-test
result_log = t.test(log10(tog$mh_pgen), log10(tog$no_mh_pgen), paired = TRUE)
result = t.test(tog$mh_pgen, tog$no_mh_pgen, paired = TRUE)


tog[, diff := log10(mh_pgen) - log10(no_mh_pgen)]

compare2 = ggplot(tog)+
    geom_hex(aes(y = diff, x = avg_mh))+
    # geom_smooth(method = 'lm', aes(x = avg_index_mh_by_pair, y = prop), size = 0.5) +
    theme_cowplot(font_family = 'Arial') + 
    ylab('log10 difference in sequence generation probability\n(MH vs no-MH model)') +
    xlab('Average MH across annotations') +
    geom_abline(intercept = 0, slope = 0, color = 'black', linewidth = 1) +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 0.5) +
    scale_fill_gradient(low = '#e8f5f1', high = '#1b9e77', name = 'count')+
    theme(text = element_text(size = 16), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 12), plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), legend.position = 'bottom', legend.justification="center", legend.key.width = unit(1, "in"), legend.key.height = unit(1, "cm"), panel.grid.major = element_line(size = 0.25))

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_annotation_ranking/pgen_diff.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = compare2, width = 8.5, height = 8.5, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)
