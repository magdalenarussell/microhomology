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

# show how probabilities differ as you increase MH
tog[, prob_diff := predicted_prob.x - predicted_prob.y]

reg_result = lm(prob_diff ~ ligation_mh, data = tog)

diffs = ggplot(tog) +
    geom_bin_2d(aes(x = as.factor(ligation_mh), y = prob_diff), binwidth = c(1, 0.005))+
    theme_cowplot(font_family = 'Arial') + 
    xlab('Number of MH nucleotides within\ntrimming and ligation configuration') +
    ylab('Difference in joint trimming and ligation probabilities\n(model with MH - model without MH)')+
    geom_abline(slope = 0, intercept = 0, color = 'black', linewidth = 0.5)+
    # geom_smooth(method = 'lm', aes(x = ligation_mh, y = prob_diff), size = 0.5) +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 0.5) +
    scale_fill_gradient(low = '#fbefe5', high = '#d95f02', trans = 'log10', name = "count") +
    theme(text = element_text(size = 8), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 6), plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), legend.position = 'bottom', legend.justification="center", legend.key.height = unit(0.3, "cm"), panel.grid.major = element_line(size = 0.25)) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_annotation_ranking/scatter_prob_diff.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = diffs, width = 8.8, height = 9.5, units = 'cm', dpi = 750, device = cairo_pdf, limitsize=FALSE)



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
index_mh_avg = tog[, mean(ligation_mh), by = .(v_gene, j_gene, index)][, mean(V1), by = .(v_gene, j_gene)]
setnames(index_mh_avg, 'V1', 'avg_index_mh_by_pair')
first_rank_subset = merge(first_rank_subset, index_mh_avg, by = c('v_gene', 'j_gene'))


# show distribution of top-ranked annotation probability differences
no_cols = c('index', 'no_annot_prob', 'no_annot_rank', 'same_top_rank', 'avg_index_mh')
no_mh_first_rank = unique(first_rank_subset[, ..no_cols])[no_annot_rank == 1]

mh_cols = c('index', 'mh_annot_prob', 'mh_annot_rank', 'same_top_rank', 'avg_index_mh')
mh_first_rank = unique(first_rank_subset[, ..mh_cols])[mh_annot_rank == 1]

rank_tog = merge(mh_first_rank, no_mh_first_rank, by = c('index', 'avg_index_mh', 'same_top_rank'))
rank_tog[, annot_diff := abs(mh_annot_prob - no_annot_prob)]

raw_prop = rank_tog[, .N, by = same_top_rank]
raw_prop[, prop:= N/sum(N)]

# get per gene prop
changed_simple = first_rank_subset[, .N, by = .(index, v_gene, j_gene, same_top_rank, avg_index_mh_by_pair, avg_pair_mh)]
s2 = changed_simple[, .N, by = .(v_gene, j_gene, same_top_rank, avg_index_mh_by_pair, avg_pair_mh)]
s2[, prop:= N/sum(N), by = .(v_gene, j_gene)]
s2$indicator = 1

# compare MH to number of changes
regression_result = lm(prop ~ avg_index_mh_by_pair, data = s2[same_top_rank == FALSE])
cor.test( s2[same_top_rank == FALSE]$prop,  s2[same_top_rank == FALSE]$avg_index_mh_by_pair, method = "pearson")

annot_changes = ggplot(s2[same_top_rank == FALSE])+
    geom_hex(aes(y = avg_index_mh_by_pair, x = prop))+
    # geom_smooth(method = 'lm', aes(x = avg_index_mh_by_pair, y = prop), size = 0.5) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('Proportion of sequences per gene pair\nwith a different top-ranked annotations\n(model with MH vs model without MH)') +
    ylab('Average MH across sequence annotation\nscenarios per gene pair') +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 0.5) +
    scale_fill_gradient(low = '#e8f5f1', high = '#1b9e77', name = 'count')+
    theme(text = element_text(size = 8), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 6), plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), legend.position = 'bottom', legend.justification="center", legend.key.height = unit(0.3, "cm"), panel.grid.major = element_line(size = 0.25)) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_annotation_ranking/per_gene_annot_changes_by_mh.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = annot_changes, width = 8.8, height = 9.5, units = 'cm', dpi = 750, device = cairo_pdf, limitsize=FALSE)

annot_change_density = ggplot(s2[same_top_rank == FALSE])+
    geom_density(aes(x = prop, y = after_stat(count)), fill = '#1b9e77', alpha = 0.8, linewidth = 0.25) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('Proportion of sequences per gene pair\nwith a different top-ranked annotations\n(model with MH vs model without MH)') +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 0.5) +
    theme(text = element_text(size = 8), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 6), plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), legend.position = 'none', panel.grid.major = element_line(size = 0.25)) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_annotation_ranking/per_gene_annot_changes_density.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = annot_change_density, width = 8.8, height = 3, units = 'cm', dpi = 750, device = cairo_pdf, limitsize=FALSE)


# explore the meaningfulness of the annotation probability differences
first_diffs = first_rank_subset[, .(abs(diff(no_annot_prob)), abs(diff(mh_annot_prob))), by = .(index, same_top_rank, avg_index_mh)]
setnames(first_diffs, 'V1', 'abs_diff_no_annot_prob')
setnames(first_diffs, 'V2', 'abs_diff_mh_annot_prob')

meaning = ggplot(first_diffs) +
    geom_hex(aes(x = abs_diff_mh_annot_prob, y = abs_diff_no_annot_prob)) +
    theme_cowplot(font_family = 'Arial') + 
    xlim(c(0, 0.825)) +
    ylim(c(0, 0.825)) +
    xlab('Absolute difference in MH model probabilities\n(for top MH vs top no-MH annotations)')+
    ylab('Absolute difference in no-MH model probabilities\n(for top MH vs top no-MH annotations)')+
    panel_border(color = 'gray60', size = 0.5) +
    background_grid(major = 'xy') +
    scale_fill_gradient(low = '#f1f0f7', high = '#7570b3', trans = 'log10', name = "count", na.value = 'white')+
    theme(text = element_text(size = 8), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 6), plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), legend.position = 'bottom', legend.justification="center", legend.key.height = unit(0.3, "cm"),panel.grid.major = element_line(size = 0.25))

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_annotation_ranking/change_magnitude_heatmap.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = meaning, width = 8.8, height = 9.5, units = 'cm', dpi = 750, device = cairo_pdf, limitsize=FALSE)

mag_density = ggplot(first_diffs)+
    geom_density(aes(x = abs_diff_mh_annot_prob, y = after_stat(count)), fill = '#7570b3', alpha = 0.8, linewidth = 0.25) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('Absolute difference in MH model probabilities\n(for top MH vs top no-MH annotations)')+
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 0.5) +
    theme(text = element_text(size = 8), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 6), plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), legend.position = 'none', panel.grid.major = element_line(size = 0.25)) 

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/analysis_annotation_ranking/change_magnitude_density.pdf')
dir.create(dirname(file_name), recursive = TRUE)

ggsave(file_name, plot = mag_density, width = 8.8, height = 3, units = 'cm', dpi = 750, device = cairo_pdf, limitsize=FALSE)
