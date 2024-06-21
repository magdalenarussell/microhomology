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
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))
NCPU <<- as.numeric(args[2])
# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 1
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 2
MODEL_TYPE <<- 'motif_two-side-base-count-beyond_average-mh_ligation-mh'
L2 <<- 'True'

source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

predictions_file = get_validation_predictions_file_path(L2, validation_annotation = 'igor_prob_experiment', sample_annotation=SAMPLE_ANNOT)

pred = fread(predictions_file)[!(new_type == '0')]

# get number of zero MH and nonzero MH annotations per index
pred[ligation_mh == 0, zeroMH_indicator := 1]
pred[ligation_mh != 0, zeroMH_indicator := 0]
pred[ligation_mh == 0, nonzeroMH_indicator := 0]
pred[ligation_mh != 0, nonzeroMH_indicator := 1]

pred[, zeroMH_annotation_count := sum(zeroMH_indicator), by = index]
pred[, nonzeroMH_annotation_count := sum(nonzeroMH_indicator), by = index]
pred[, average_annotation_MH := mean(ligation_mh), by = index]
pred[, total_annotations := .N, by = index]

# normalize probabilities
igor_probs = pred[ligation_mh == 0]
igor_probs[, igor_sum := sum(norm_default_igor_prob, na.rm = TRUE), by = .(v_gene, j_gene)]
igor_probs[, norm_igor_prob := norm_default_igor_prob/igor_sum]

# remove indices which do not have at least one zeroMH annotation (these are
# all edge cases near trimming domain restriction)
pred_subset = pred[!(zeroMH_annotation_count == 0)]
igor_probs_subset = igor_probs[!(zeroMH_annotation_count == 0)]

# aggregate probabilities by index
cols = c('v_gene', 'j_gene', 'seq_index', 'zeroMH_annotation_count', 'nonzeroMH_annotation_count', 'average_annotation_MH', 'total_annotations')
pred_agg = pred_subset[, sum(predicted_prob, na.rm = TRUE), by = cols]
setnames(pred_agg, 'V1', 'summed_norm_predicted_prob')

# get per-zero-annotation estimated probability by index
pred_agg[, perseq_norm_predicted_prob := summed_norm_predicted_prob/zeroMH_annotation_count]

# combine
tog_subset = merge(igor_probs_subset, pred_agg, by = cols)

# bin total annotations and nonzeroMH annotation count
tog_subset[total_annotations >= 12, binned_total_annotations := '12 or more']
tog_subset[total_annotations < 12, binned_total_annotations := as.character(total_annotations)]
tog_subset[nonzeroMH_annotation_count < 8, binned_nonzeroMH_annotation_count := as.character(nonzeroMH_annotation_count)]
tog_subset[nonzeroMH_annotation_count >= 8, binned_nonzeroMH_annotation_count := '8 or more']

# get gene probabilities from IGoR (filter out low prob genes)
j_usage = fread(paste0(MOD_OUTPUT_PATH, '/meta_data/igor_alpha_jchoice_params.tsv'))
v_usage = fread(paste0(MOD_OUTPUT_PATH, '/meta_data/igor_alpha_vchoice_params.tsv'))
j_usage$indicator = 1
v_usage$indicator = 1

joint_usage = merge(v_usage, j_usage, by = 'indicator', allow.cartesian = TRUE)
joint_usage[, joint_gene_usage := v_gene_prob * j_gene_prob]

## merge and restrict to genes that have high usage in IGoR
tog_subset = merge(tog_subset, joint_usage, by = c('v_gene', 'j_gene'))
tog_subset2 = tog_subset[joint_gene_usage > 1e-4]

plot = ggplot(tog_subset2) +
    geom_point(aes(x = norm_igor_prob, y = perseq_norm_predicted_prob), size = 2, alpha = 0.1) +
    geom_abline(intercept = 0, slope = 1, color = 'blue', linewidth = 1) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('Normalized IGoR probability (default TRA model)') +
    ylab('MH model probability') +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_igor_comparison/default_TRA_igor_comparison.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = plot, width = 10, height = 8, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)

plot_split = ggplot(tog_subset2) +
    geom_point(aes(x = norm_igor_prob, y = perseq_norm_predicted_prob, color = zeroMH_annotation_count), size = 2, alpha = 0.3) +
    facet_wrap(~zeroMH_annotation_count, nrow = 3, ncol = 3) +
    geom_abline(intercept = 0, slope = 1, color = 'blue', linewidth = 1) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('Normalized IGoR probability (default TRA model)') +
    ylab('MH model probability') +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_igor_comparison/default_TRA_igor_comparison_zeroMH.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = plot_split, width = 20, height = 17, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)


plot_split2 = ggplot(tog_subset2) +
    geom_point(aes(x = norm_igor_prob, y = perseq_norm_predicted_prob, color = binned_nonzeroMH_annotation_count), size = 2, alpha = 0.3) +
    facet_wrap(~binned_nonzeroMH_annotation_count, nrow = 3, ncol = 3) +
    geom_abline(intercept = 0, slope = 1, color = 'blue', linewidth = 1) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('Normalized IGoR probability (default TRA model)') +
    ylab('MH model probability') +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_igor_comparison/default_TRA_igor_comparison_nonzeroMH.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = plot_split2, width = 20, height = 17, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)


plot2 = ggplot(tog_subset2) +
    geom_point(aes(x = log(norm_igor_prob), y = log(perseq_norm_predicted_prob)), size = 2, alpha = 0.1) +
    geom_abline(intercept = 0, slope = 1, color = 'blue', linewidth = 1) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('(log) Normalized IGoR probability (default TRA model)') +
    ylab('(log) MH model probability') +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_igor_comparison/default_TRA_igor_comparison_log.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = plot2, width = 10, height = 8, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)

plot2_split = ggplot(tog_subset2) +
    geom_point(aes(x = log(norm_igor_prob), y = log(perseq_norm_predicted_prob), color = binned_nonzeroMH_annotation_count), size = 2, alpha = 0.3) +
    facet_wrap(~binned_nonzeroMH_annotation_count, nrow = 3, ncol = 3) +
    geom_abline(intercept = 0, slope = 1, color = 'blue', linewidth = 1) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('Normalized IGoR probability (default TRA model)') +
    ylab('MH model probability') +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_igor_comparison/default_TRA_igor_comparison_nonzeroMH_log.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = plot2_split, width = 20, height = 17, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)

plot2_split2 = ggplot(tog_subset2) +
    geom_point(aes(x = log(norm_igor_prob), y = log(perseq_norm_predicted_prob), color = zeroMH_annotation_count), size = 2, alpha = 0.3) +
    facet_wrap(~zeroMH_annotation_count, nrow = 3, ncol = 3) +
    geom_abline(intercept = 0, slope = 1, color = 'blue', linewidth = 1) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('Normalized IGoR probability (default TRA model)') +
    ylab('MH model probability') +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_igor_comparison/default_TRA_igor_comparison_zeroMH_log.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = plot2_split2, width = 20, height = 17, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)

plot3_split = ggplot(tog_subset2) +
    geom_point(aes(x = norm_igor_prob, y = perseq_norm_predicted_prob, color = average_annotation_MH), size = 2, alpha = 0.3) +
    facet_wrap(~binned_total_annotations, nrow = 4, ncol = 3) +
    geom_abline(intercept = 0, slope = 1, color = 'blue', linewidth = 1) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('Normalized IGoR probability (default TRA model)') +
    ylab('MH model probability') +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_igor_comparison/default_TRA_igor_comparison_avgMH.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = plot3_split, width = 20, height = 22, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)

plot3_split2 = ggplot(tog_subset2) +
    geom_point(aes(x = log(norm_igor_prob), y = log(perseq_norm_predicted_prob), color = average_annotation_MH), size = 2, alpha = 0.3) +
    facet_wrap(~binned_total_annotations, nrow = 4, ncol = 3) +
    geom_abline(intercept = 0, slope = 1, color = 'blue', linewidth = 1) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('Normalized IGoR probability (default TRA model)') +
    ylab('MH model probability') +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_igor_comparison/default_TRA_igor_comparison_avgMH_log.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = plot3_split2, width = 20, height = 22, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)

tog_subset2[, log_diff := log(perseq_norm_predicted_prob) - log(norm_igor_prob)]

plot_trend = ggplot(tog_subset2) +
    geom_hex(aes(x = average_annotation_MH, y = log_diff), binwidth = c(0.15, 2)) +
    geom_abline(intercept = 0, slope = 0, color = 'black', linewidth = 1) +
    geom_smooth(aes(x = average_annotation_MH, y = log_diff), method = 'lm', se = TRUE, linewidth = 1) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('Average MH per annotation') +
    ylab('log(MH model probability) - log(IGoR probability)') +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    scale_fill_viridis_c(option = "C", trans = "log10", name = 'Log-scaled counts') +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_igor_comparison/default_TRA_igor_comparison_trend_avgMH.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = plot_trend, width = 10, height = 8, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)

plot_trend2 = ggplot(tog_subset2[zeroMH_annotation_count == 1]) +
    geom_point(aes(x = average_annotation_MH, y = log_diff), size = 2, alpha = 0.1) +
    geom_abline(intercept = 0, slope = 0, color = 'orange', linewidth = 1) +
    geom_smooth(aes(x = average_annotation_MH, y = log_diff), method = 'lm', se = TRUE, linewidth = 1) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('Average MH per annotation') +
    ylab('log(MH model probability) - log(IGoR probability)') +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_igor_comparison/default_TRA_igor_comparison_trend_avgMH_only_one_zeroMH_annotation.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = plot_trend2, width = 10, height = 8, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)

plot_trend3 = ggplot(tog_subset2[v_gene_j_gene == 'TRAV19*01_TRAJ45*01']) +
    geom_point(aes(x = average_annotation_MH, y = log_diff), size = 2, alpha = 0.1) +
    geom_abline(intercept = 0, slope = 0, color = 'orange', linewidth = 1) +
    geom_smooth(aes(x = average_annotation_MH, y = log_diff), method = 'lm', se = TRUE, linewidth = 1) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('Average MH per annotation') +
    ylab('log(MH model probability) - log(IGoR probability)') +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_igor_comparison/default_TRA_igor_comparison_trend_avgMH_TRAV19_TRAJ45.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = plot_trend3, width = 10, height = 8, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)

plot_trend4 = ggplot(tog_subset2[v_trim >= 0 & v_trim < 9][j_trim >= 0 & j_trim <9]) +
    geom_point(aes(x = average_annotation_MH, y = log_diff), size = 2, alpha = 0.1) +
    geom_abline(intercept = 0, slope = 0, color = 'orange', linewidth = 1) +
    geom_smooth(aes(x = average_annotation_MH, y = log_diff), method = 'lm', se = TRUE, linewidth = 1) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('Average MH per annotation') +
    ylab('log(MH model probability) - log(IGoR probability)') +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_igor_comparison/default_TRA_igor_comparison_trend_avgMH_trims0_to_8.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = plot_trend4, width = 10, height = 8, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)
