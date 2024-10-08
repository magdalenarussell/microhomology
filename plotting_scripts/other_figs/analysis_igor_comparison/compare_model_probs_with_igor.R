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

# get igor probabilities
# get IGoR Vtrim and Jtrim probs
jt = fread(paste0(MOD_OUTPUT_PATH, '/igor_alpha/nonproductive_v-j_trim_ligation-mh/unbounded_motif_trims_bounded_-2_14/1_2_motif_two-side-base-count-beyond_average-mh_ligation-mh/igor_prob_experiment/igor_alpha_j_trim_params.tsv'))
vt = fread(paste0(MOD_OUTPUT_PATH, '/igor_alpha/nonproductive_v-j_trim_ligation-mh/unbounded_motif_trims_bounded_-2_14/1_2_motif_two-side-base-count-beyond_average-mh_ligation-mh/igor_prob_experiment/igor_alpha_v_trim_params.tsv'))
jt$indicator = 1
vt$indicator = 1

igor_probs = merge(vt, jt, by = "indicator", allow.cartesian = TRUE)
igor_probs[, igor_joint_prob := v_trim_prob * j_trim_prob]
igor_probs[, igor_sum := sum(igor_joint_prob, na.rm = TRUE), by = .(v_gene, j_gene)]
igor_probs[, norm_igor_joint_prob := igor_joint_prob/igor_sum]

# join everything together
cols = c('index', 'v_gene', 'j_gene', 'v_trim', 'j_trim', 'ligation_mh')
pcols = c('index', 'v_gene', 'j_gene', 'v_trim', 'j_trim', 'ligation_mh', 'processed_sequence', 'predicted_prob')
icols = c('v_gene', 'j_gene', 'v_trim', 'j_trim', 'ligation_mh', 'norm_igor_joint_prob')

igor_probs$ligation_mh = 0
igor_tog = merge(igor_probs[, ..icols], unique(pred[ligation_mh == 0][, ..cols]), by = c('v_gene', 'j_gene', 'v_trim', 'j_trim', 'ligation_mh'))

# .x vars are mh mapping and .y vars are igor mapping
tog = merge(pred[, ..pcols], igor_tog, by = c('index', 'v_gene', 'j_gene'), allow.cartesian = TRUE)
tog = tog[!((ligation_mh.x == 0 & ligation_mh.y == 0) & (v_trim.x != v_trim.y | j_trim.x != j_trim.y))]

# get the number of mh annotations that map to each igor annotation
tog[, mh_to_igor_annot_count := .N, by = .(v_trim.y, j_trim.y, ligation_mh.y, index, v_gene, j_gene)]

# get the number of igor annotations that map to each mh annotation
tog[, igor_to_mh_annot_count := .N, by = .(v_trim.x, j_trim.x, ligation_mh.x, index, v_gene, j_gene)]

# Adjust probabilities by weighting them by mapping count
tog[, adjusted_predicted_prob := predicted_prob/igor_to_mh_annot_count]
tog[, adjusted_norm_igor_joint_prob := norm_igor_joint_prob/mh_to_igor_annot_count]

# Aggregate MH annotation probabilities from IGOR probabilities
mh_annot = tog[, sum(adjusted_norm_igor_joint_prob), by = .(index, v_gene, j_gene, v_trim.x, j_trim.x, ligation_mh.x, igor_to_mh_annot_count)]
mh_annot[, sum_converted_igor_prob := sum(V1), by = .(v_gene, j_gene)]
mh_annot[, norm_converted_mh_annot_igor_prob := V1/sum_converted_igor_prob]

# merge with mh model predictions
mh_tog = merge(pred, mh_annot, by.x = c('index', 'v_gene', 'j_gene', 'v_trim', 'j_trim', 'ligation_mh'), by.y = c('index', 'v_gene', 'j_gene', 'v_trim.x', 'j_trim.x', 'ligation_mh.x'))

# Aggregate IGOR annotation probabilities from MH probabilities
tog[ligation_mh.x == 0, nonzeroMH := 0]
tog[ligation_mh.x > 0, nonzeroMH := 1]
tog[, sum_mh_configs := sum(nonzeroMH), by = .(index, v_gene, j_gene, v_trim.y, j_trim.y)]
tog[, sum_mh := sum(ligation_mh.x), by = .(index, v_gene, j_gene, v_trim.y, j_trim.y)]
igor_annot = tog[, sum(adjusted_predicted_prob), by = .(index, v_gene, j_gene, v_trim.y, j_trim.y, ligation_mh.y, mh_to_igor_annot_count, sum_mh, sum_mh_configs)]
igor_annot[, avg_mh := sum_mh/mh_to_igor_annot_count]
igor_annot[, mh_config_prop := sum_mh_configs/mh_to_igor_annot_count]
igor_annot[, sum_converted_mh_prob := sum(V1), by = .(v_gene, j_gene)]
igor_annot[, norm_converted_igor_annot_mh_prob := V1/sum_converted_mh_prob]

## merge with igor model predictions
igor_tog2 = merge(igor_tog, igor_annot, by.x = c('index', 'v_gene', 'j_gene', 'v_trim', 'j_trim', 'ligation_mh'), by.y = c('index', 'v_gene', 'j_gene', 'v_trim.y', 'j_trim.y', 'ligation_mh.y'))

# subset by gene usage
# get gene probabilities from IGoR (filter out low prob genes)
j_usage = fread(paste0(MOD_OUTPUT_PATH, '/igor_alpha/nonproductive_v-j_trim_ligation-mh/unbounded_motif_trims_bounded_-2_14/1_2_motif_two-side-base-count-beyond_average-mh_ligation-mh/igor_prob_experiment/igor_alpha_jchoice_params.tsv'))
v_usage = fread(paste0(MOD_OUTPUT_PATH, '/igor_alpha/nonproductive_v-j_trim_ligation-mh/unbounded_motif_trims_bounded_-2_14/1_2_motif_two-side-base-count-beyond_average-mh_ligation-mh/igor_prob_experiment/igor_alpha_vchoice_params.tsv'))
j_usage$indicator = 1
v_usage$indicator = 1

joint_usage = merge(v_usage, j_usage, by = 'indicator', allow.cartesian = TRUE)
joint_usage[, joint_gene_usage := v_gene_prob * j_gene_prob]

## merge and restrict to genes that have high usage in IGoR
igor_tog2_merge = merge(igor_tog2, joint_usage, by = c('v_gene', 'j_gene'))
igor_tog_subset = igor_tog2_merge[joint_gene_usage > 1e-4]

mh_tog_merge = merge(mh_tog, joint_usage, by = c('v_gene', 'j_gene'))
mh_tog_subset = mh_tog_merge[joint_gene_usage > 1e-4]

# make sure all probabilities are normalized
mh_tog_subset[, predicted_prob2 := predicted_prob/sum(predicted_prob), by = .(v_gene, j_gene)]
igor_tog_subset[, norm_igor_joint_prob2 := norm_igor_joint_prob/sum(norm_igor_joint_prob), by = .(v_gene, j_gene)]

# get diff
mh_tog_subset[, diff := predicted_prob2 - norm_converted_mh_annot_igor_prob]
igor_tog_subset[, diff := norm_converted_igor_annot_mh_prob - norm_igor_joint_prob2]

# Perform correlation analysis
mh_cor = cor(mh_tog_subset$diff, mh_tog_subset$ligation_mh)
igor_cor = cor(igor_tog_subset$diff, igor_tog_subset$avg_mh)
mh_reg = lm(diff ~ ligation_mh + igor_to_mh_annot_count, data = mh_tog_subset)
igor_reg = lm(diff ~ avg_mh + mh_to_igor_annot_count, data = igor_tog_subset)

# compare MH annotation probabilities
plot_trend = ggplot(mh_tog_subset) +
    geom_hex(aes(x = predicted_prob2, y = norm_converted_mh_annot_igor_prob)) +
    geom_abline(intercept = 0, slope = 1, color = 'black', linewidth = 1) +
    theme_cowplot(font_family = 'Arial') + 
    ylab("Aggregated IGoR probability\n(MH annotation space)") +
    xlab("MH model probability") +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    scale_fill_gradient(low = '#efedf5', high = '#3f007d', trans = 'log10', name = expression(log[10]("count"))) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_igor_comparison/TRA_igor_comparison_mh_annotation.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = plot_trend, width = 10, height = 7, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)

plot_trend = ggplot(mh_tog_subset) +
    geom_hex(aes(x = ligation_mh, y = diff)) +
    geom_abline(intercept = 0, slope = 0, color = 'black', linewidth = 1) +
    geom_smooth(aes(x = ligation_mh, y = diff), method = 'lm', se = TRUE, linewidth = 1) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('MH within trimming/ligation annotation') +
    ylab("MH model probability - aggregated IGoR probability\n(MH annotation space)")+
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    scale_fill_gradient(low = '#efedf5', high = '#3f007d', trans = 'log10', name = expression(log[10]("count"))) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_igor_comparison/TRA_igor_comparison_mh_annotation_ligMH.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = plot_trend, width = 10, height = 7, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)



# compare IGoR annotation probabilities
plot_trend = ggplot(igor_tog_subset) +
    geom_hex(aes(x = norm_converted_igor_annot_mh_prob, y = norm_igor_joint_prob2)) +
    geom_abline(intercept = 0, slope = 1, color = 'black', linewidth = 1) +
    theme_cowplot(font_family = 'Arial') + 
    ylab("IGoR probability")+
    xlab("Aggregated MH model probability\n(IGoR annotation space)") +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    scale_fill_gradient(low = '#efedf5', high = '#3f007d', trans = 'log10', name = expression(log[10]("count"))) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_igor_comparison/TRA_igor_comparison_igor_annotation.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = plot_trend, width = 10, height = 7, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)

plot_trend = ggplot(igor_tog_subset) +
    geom_hex(aes(x = avg_mh, y = diff)) +
    geom_abline(intercept = 0, slope = 0, color = 'black', linewidth = 1) +
    geom_smooth(aes(x = avg_mh, y = diff), method = 'lm', se = TRUE, linewidth = 1) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('Average MH within trimming/ligation annotations') +
    ylab("Aggregated MH model probability - IGoR probability\n(IGoR annotation space)")+
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    scale_fill_gradient(low = '#efedf5', high = '#3f007d', trans = 'log10', name = expression(log[10]("count"))) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_igor_comparison/TRA_igor_comparison_igor_annotation_ligMH.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = plot_trend, width = 10, height = 7, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)




