source('config/config.R')

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

L2 <<- 'False'

source(paste0(MOD_PROJECT_PATH,'/analysis_scripts/ligation-mh_signal_simulator_scripts/ligation-mh_simulator_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

# get all possible configurations
configs = read_frames_data()

configs[ligation_mh == 0, mh_type := 'zero']
configs[ligation_mh > 0, mh_type := 'nonzero']

configs[, all_options := .N, by = .(v_gene, j_gene, processed_sequence)]
configs[, mh_options := .N, by = .(v_gene, j_gene, processed_sequence, mh_type)]

cols = c('v_gene', 'j_gene', 'processed_sequence', 'processed_nt_change', 'mh_type', 'all_options', 'mh_options')
subset = unique(configs[, ..cols])

subset_wide = subset %>% pivot_wider(names_from = 'mh_type', values_from = 'mh_options') %>% as.data.table()
subset_wide[is.na(zero), zero := 0]
subset_wide[is.na(nonzero), nonzero := 0]

subset_wide_condensed = subset_wide[, .N, by = .(zero, nonzero)]
subset_wide_condensed = subset_wide_condensed[!(zero == 0 & nonzero > 0)]

plot = ggplot(subset_wide_condensed, aes(x = zero, y = nonzero, fill = N)) + 
       geom_tile() +
       xlab("Zero MH annotation count") +
       ylab("Nonzero MH annotation count") +
       theme_cowplot(font_family = 'Arial') + 
       background_grid(major = 'xy') + 
       scale_fill_distiller(palette = 'YlOrRd', trans = 'log', breaks = c(1, 10, 100, 1000, 100000), labels = c("1", "10", "100", "1,000", "100,000")) +
       panel_border(color = 'gray60', size = 1.5) +
       theme(text = element_text(size = 40), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 30), legend.key.size = unit(3, "cm")) +
       labs('Sequence count')

# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/ligation_mh_choices_germline_tra/annotation_choices_scatter.pdf')

ggsave(file_name, plot = plot, width = 18, height = 16, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE) 

subset_wide_type = subset_wide[!(zero == 0 & nonzero > 0), median(processed_nt_change), by = .(zero, nonzero)] 

plot = ggplot(subset_wide_type, aes(x = zero, y = nonzero, fill = V1)) + 
       geom_tile() +
       xlab("Zero MH annotation count") +
       ylab("Nonzero MH annotation count") +
       theme_cowplot(font_family = 'Arial') + 
       background_grid(major = 'xy') + 
       scale_fill_distiller(palette = 'YlOrRd') +
       panel_border(color = 'gray60', size = 1.5) +
       theme(text = element_text(size = 40), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 30), legend.key.size = unit(3, "cm")) +
       labs(fill = 'Number of processed\nnucleotides')

# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/ligation_mh_choices_germline_tra/annotation_nt_processing.pdf')

ggsave(file_name, plot = plot, width = 18, height = 16, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE) 


plot = ggplot(subset) + 
       geom_histogram(aes(x = all_options), fill = 'gray', color = 'gray', binwidth = 2.5, alpha = 0.8) +
       # geom_histogram(aes(x = mh_options, fill = mh_type, color = mh_type), alpha = 0.5) +
       xlab("Annotation count") +
       theme_cowplot(font_family = 'Arial') + 
       background_grid(major = 'xy') + 
       panel_border(color = 'gray60', size = 1.5) +
       theme(text = element_text(size = 40), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 30), legend.key.size = unit(3, "cm")) 
       # labs('Amount of MH in annotation')

# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/ligation_mh_choices_germline_tra/annotation_choices_density.pdf')

ggsave(file_name, plot = plot, width = 16, height = 5, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE) 

trimming_probs = get_igor_trimming_probs()

configs = merge(configs, trimming_probs$v_trim[, -c('prob_sum')], by = c('v_gene', 'v_trim'))
configs = merge(configs, trimming_probs$j_trim[, -c('prob_sum')], by = c('j_gene', 'j_trim'))

configs[, joint_trimming_prob := v_trim_prob * j_trim_prob]
configs[joint_trimming_prob == 0, joint_trimming_prob := 1e-20]

configs[, joint_trimming_prob_norm := joint_trimming_prob/sum(joint_trimming_prob), by = .(v_gene, j_gene, processed_sequence)]
mh_trim = configs[, median(joint_trimming_prob_norm), by = .(v_gene,  j_gene, processed_sequence, mh_type, processed_nt_change)]

mh_trim_wide = mh_trim %>% pivot_wider(names_from = 'mh_type', values_from = 'V1') %>% as.data.table()

plot = ggplot(mh_trim_wide, aes(x = zero, y = nonzero, color = processed_nt_change)) + 
       geom_point(size = 4, alpha = 0.5) +
       xlab("Zero MH annotation\nmedian joint trimming probability") +
       ylab("Nonzero MH annotation count\nmedian joint trimming probability") +
       theme_cowplot(font_family = 'Arial') + 
       background_grid(major = 'xy') + 
       scale_color_distiller(palette = 'YlOrRd') +
       panel_border(color = 'gray60', size = 1.5) +
       theme(text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 22), legend.key.size = unit(3, "cm")) +
       labs(color = 'Number of processed\nnucleotides') 

# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/ligation_mh_choices_germline_tra/annotation_trim_probs.pdf')

ggsave(file_name, plot = plot, width = 18, height = 14, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE) 

plot = ggplot(mh_trim_wide, aes(x = log(zero), y = log(nonzero), color = processed_nt_change)) + 
       geom_point(size = 4, alpha = 0.5) +
       xlab("Zero MH annotation count\nmedian joint trimming probability (log)") +
       ylab("Nonzero MH annotation count\nmedian joint trimming probability (log)") +
       theme_cowplot(font_family = 'Arial') + 
       background_grid(major = 'xy') + 
       scale_color_distiller(palette = 'YlOrRd') +
       panel_border(color = 'gray60', size = 1.5) +
       theme(text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 22), legend.key.size = unit(3, "cm")) +
       labs(color = 'Number of processed\nnucleotides') 

# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/ligation_mh_choices_germline_tra/annotation_trim_probs_log.pdf')

ggsave(file_name, plot = plot, width = 18, height = 14, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE) 


