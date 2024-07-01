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

# re-normalize
pred_subset[, pred_sum := sum(predicted_prob), by = .(v_gene, j_gene)]
igor_probs_subset[, igor_sum := sum(norm_igor_prob), by = .(v_gene, j_gene)]

pred_subset[, predicted_prob := predicted_prob/pred_sum]
igor_probs_subset[, norm_igor_prob := norm_igor_prob/igor_sum]

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
j_usage = fread(paste0(MOD_OUTPUT_PATH, '/igor_alpha/nonproductive_v-j_trim_ligation-mh/unbounded_motif_trims_bounded_-2_14/1_2_motif_two-side-base-count-beyond_average-mh_ligation-mh/igor_prob_experiment/igor_alpha_jchoice_params.tsv'))
v_usage = fread(paste0(MOD_OUTPUT_PATH, '/igor_alpha/nonproductive_v-j_trim_ligation-mh/unbounded_motif_trims_bounded_-2_14/1_2_motif_two-side-base-count-beyond_average-mh_ligation-mh/igor_prob_experiment/igor_alpha_vchoice_params.tsv'))
j_usage$indicator = 1
v_usage$indicator = 1

joint_usage = merge(v_usage, j_usage, by = 'indicator', allow.cartesian = TRUE)
joint_usage[, joint_gene_usage := v_gene_prob * j_gene_prob]

## merge and restrict to genes that have high usage in IGoR
tog_subset = merge(tog_subset, joint_usage, by = c('v_gene', 'j_gene'))
tog_subset2 = tog_subset[joint_gene_usage > 1e-4]

tog_subset2[, log_diff := log10(perseq_norm_predicted_prob) - log10(norm_igor_prob)]
tog_subset2[, diff := perseq_norm_predicted_prob - norm_igor_prob]


plot_trend = ggplot(tog_subset2) +
    geom_hex(aes(x = average_annotation_MH, y = log_diff), binwidth = c(0.14, 5)) +
    # geom_hex(aes(x = average_annotation_MH, y = log_diff)) +
    geom_abline(intercept = 0, slope = 0, color = 'black', linewidth = 1) +
    geom_smooth(aes(x = average_annotation_MH, y = log_diff), method = 'lm', se = TRUE, linewidth = 1) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('Average MH per annotation') +
    ylab(expression(log[10]("MH model probability") - log[10]("IGoR probability")))+
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    scale_fill_gradient(low = '#efedf5', high = '#3f007d', trans = 'log10', name = expression(log[10]("count"))) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_igor_comparison/TRA_igor_comparison_trend_avgMH.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = plot_trend, width = 10, height = 7, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)

plot_trend = ggplot(tog_subset2) +
    geom_hex(aes(x = average_annotation_MH, y = diff)) +
    geom_abline(intercept = 0, slope = 0, color = 'black', linewidth = 1) +
    geom_smooth(aes(x = average_annotation_MH, y = diff), method = 'lm', se = TRUE, linewidth = 1) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('Average MH per annotation') +
    ylab("MH model probability - IGoR probability")+
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    scale_fill_gradient(low = '#efedf5', high = '#3f007d', trans = 'log10', name = expression(log[10]("count"))) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_igor_comparison/TRA_igor_comparison_trend_avgMH_nonlog.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = plot_trend, width = 10, height = 7, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)
