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

# join everything together .x vars are mh model and .y vars are noMH model
tog = merge(pred, no_pred, by = c('index', 'v_gene', 'j_gene', 'v_trim', 'j_trim', 'ligation_mh'))

# convert probabilities to be p_annot
tog[, mh_model_prob_sum := sum(predicted_prob.x), by = index]
tog[, no_model_prob_sum := sum(predicted_prob.y), by = index]
tog[, mh_annot_prob := predicted_prob.x/mh_model_prob_sum]
tog[, no_annot_prob := predicted_prob.y/no_model_prob_sum]

# get annotation ranking
tog[, mh_annot_rank := rank(-mh_annot_prob), by = index]
tog[, no_annot_rank := rank(-no_annot_prob), by = index]
tog[, total_annot_count := .N, by = index]

# Explore how often the first rank changes
first_rank = tog[mh_annot_rank == 1 | no_annot_rank == 1]
first_rank_subset = first_rank[total_annot_count > 1]
first_rank_subset[mh_annot_rank == no_annot_rank, same_top_rank := TRUE]
first_rank_subset[mh_annot_rank != no_annot_rank, same_top_rank := FALSE]

# get MH
first_rank_subset[, avg_pair_mh := mean(ligation_mh), by = .(v_gene, j_gene)]
first_rank_subset[, avg_index_mh := mean(ligation_mh), by = .(v_gene, j_gene, index)]

# get per gene prop
changed_simple = first_rank_subset[, .N, by = .(index, v_gene, j_gene, same_top_rank, avg_pair_mh)]
s2 = changed_simple[, .N, by = .(v_gene, j_gene, same_top_rank, avg_pair_mh)]
s2[, prop:= N/sum(N), by = .(v_gene, j_gene)]
s2$indicator = 1

# plot side by side
s2$model = "MH model versus\nmodel without MH terms"

plot = ggplot(s2[same_top_rank == FALSE])+
    geom_boxplot(aes(x = model, y = prop, fill = model))+
    theme_cowplot(font_family = 'Arial') + 
    xlab('') +
    ylab('Proportion of sequences with different\ntop ranked MH annotation (per gene pair)\n') +
    background_grid(major = 'xy') + 
    scale_fill_brewer(palette = 'Dark2') +
    panel_border(color = 'gray60', size = 1.5) +
    theme(text = element_text(size = 14), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), legend.position = 'none') 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_annotation_ranking/boxplot_MH_noMH.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = plot, width = 5, height = 7, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)

first_diffs = first_rank_subset[, .(abs(diff(no_annot_prob)), abs(diff(mh_annot_prob))), by = .(index, same_top_rank, avg_index_mh)]
setnames(first_diffs, 'V1', 'abs_diff_no_annot_prob')
setnames(first_diffs, 'V2', 'abs_diff_mh_annot_prob')


scatter = ggplot(first_diffs) +
    geom_hex(aes(x = abs_diff_mh_annot_prob, y = abs_diff_no_annot_prob))+
    theme_cowplot(font_family = 'Arial') + 
    xlab('Absolute difference in MH model probabilities\n(between top MH model annotation and top no-MH model annotation)') +
    ylab('Absolute difference in no-MH model probabilities\n(between top MH model annotation and top no-MH model annotation)') +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    scale_fill_gradient(low = '#fbefe5', high = '#d95f02', trans = 'log10', name = expression(log[10]("count"))) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_annotation_ranking/scatter_MH_noMH.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = scatter, width = 12, height = 10, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)

# compare MH to number of changes
scatter2 = ggplot(s2[same_top_rank == FALSE])+
    geom_hex(aes(x = avg_pair_mh, y = prop))+
    geom_smooth(method = 'lm', aes(x = avg_pair_mh, y = prop)) +
    theme_cowplot(font_family = 'Arial') + 
    ylab('Proportion of sequences with different\ntop ranked MH annotation (per gene pair)\n') +
    xlab('Average MH within gene pair trimming/ligation configurations') +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    scale_fill_gradient(low = '#fbefe5', high = '#d95f02', name = 'count')+
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_annotation_ranking/scatter_MH_noMH_avgMH_effect.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = scatter2, width = 12, height = 10, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)


# compare MH to meaningfulness of changes
scatter = ggplot(first_diffs) +
    geom_hex(aes(x = avg_index_mh, y = abs_diff_mh_annot_prob))+
    geom_smooth(method = 'lm', aes(x = avg_index_mh, y = abs_diff_mh_annot_prob)) +
    xlab('Average MH across sequence annotations') +
    theme_cowplot(font_family = 'Arial') + 
    ylab('Absolute difference in MH model probabilities\n(between top MH model annotation and top no-MH model annotation)') +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    scale_fill_gradient(low = '#fbefe5', high = '#d95f02', trans = 'log10', name = expression(log[10]("count"))) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_annotation_ranking/scatter_MH_noMH_avgMH_index_effect.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = scatter, width = 12, height = 10, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)
