source('config/config.R')

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(Biostrings)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(mclogit)
library(matrixcalc)
library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)

PARAM_root <<- c('both', 'nonproductive')

NCPU <<- as.numeric(2)

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(1)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(2)

MODEL_TYPE <<- 'motif_two-side-base-count-beyond_mh-config-count_ligation-mh'

L2 <<- 'False'

LIGATION_PARAMS <<- c(0, 1)
INT_MH_PARAMS <<- c(0, 1)
TRIMMING_PROB_MODEL <<- c('motif_two-side-base-count-beyond')
SAMPLE = c(TRUE, FALSE)

results = data.table()

for (param_r in PARAM_root){
    PARAM_GROUP <<- paste0(param_r, '_v-j_trim_ligation-mh')
    source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))
    for (s in SAMPLE){
        for (trim_model in TRIMMING_PROB_MODEL){
            for (trim_mh in INT_MH_PARAMS){
                for (lig_mh in LIGATION_PARAMS){
                    ANNOTATION_TYPE <<- 'igor_mh_sim_alpha'
                    old_annotation = ANNOTATION_TYPE

                    source(paste0(MOD_PROJECT_PATH, '/config/file_paths.R'))
                    source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

                    ANNOTATION_TYPE <<- paste0(old_annotation, '_from_', trim_model, '_MHprob', lig_mh, '_trimMHprob', trim_mh)
                    # Read in model coefficient data 
                    coef_path = get_model_coef_file_path(L2, sample_annotation=s)
                    if (!file.exists(coef_path)){
                        assign(paste0('together', trim_mh, lig_mh), NULL)
                        next
                    }
                    coefs = fread(coef_path)
                    temp = coefs[coefficient %like% 'mh']
                    temp$sequence_type = param_r
                    temp$trimming_model = trim_model
                    temp$lig_MH_param = lig_mh
                    temp$trim_MH_param = trim_mh
                    if (!(coef_path %like% 'all_annotations')){
                        temp$model_type = "Coefficients from model trained without EM\n(Known annotations)"
                    } else {
                        temp$model_type = "Coefficients from model trained with EM\n(Unknown annotations)"
                    }

                    results = rbind(results, temp)
                }
            }
        }
    }
}

results[lig_MH_param == 0 & trim_MH_param == 0, simulator_type := "No MH effect for trimming or ligation"]
results[lig_MH_param == 1 & trim_MH_param == 0, simulator_type := "MH effect for ligation, but not trimming"]
results[lig_MH_param == 0 & trim_MH_param == 1, simulator_type := "MH effect for trimming, but not ligation"]
results[lig_MH_param == 1 & trim_MH_param == 1, simulator_type := "MH effect for both trimming and ligation"]

results$log_value = results$value/log(10)
cols = c('coefficient', 'log_value', 'sequence_type', 'trimming_model', 'lig_MH_param', 'trim_MH_param', 'model_type', 'simulator_type')
results_wide = results[, ..cols][coefficient %like% 'ligation'] %>% pivot_wider(names_from = model_type, values_from = log_value) %>% as.data.table()

# temp
if ('Coefficients from model trained with EM\n(Unknown annotations)' %in% colnames(results_wide)) {
    results_wide[is.na(get('Coefficients from model trained with EM\n(Unknown annotations)')), paste('Coefficients from model trained with EM\n(Unknown annotations)') := Inf]
} else {
    results_wide[['Coefficients from model trained with EM\n(Unknown annotations)']] = Inf
}

colors <- brewer.pal(4, "Set2")
plot = ggplot(results_wide, aes(x = get('Coefficients from model trained without EM\n(Known annotations)'), y = get('Coefficients from model trained with EM\n(Unknown annotations)'), color = factor(simulator_type), shape = sequence_type)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray60", linewidth=1) +
  geom_point(size = 7) +
  scale_color_manual(values = colors) +
  labs(color = "MH simulation type", shape = "Sequence/Productivity type", x = 'Coefficients from model trained without EM\n(Known annotations)', y = 'Coefficients from model trained with EM\n(Unknown annotations)') +
  theme_cowplot(font_family = 'Arial') + 
  background_grid(major = 'xy') + 
  panel_border(color = 'gray60', size = 1.1) +
  theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 18), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), panel.spacing = unit(3, "lines"), legend.key.height = unit(2, 'cm'))  

# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_simulator_experiments/compare_EM_noEM_simulator/compare_lig_motif.pdf')

ggsave(file_name, plot = plot, width = 14, height = 8, units = 'in', dpi = 750, device = cairo_pdf)

cols = c(cols, 'position')
results_wide = results[, ..cols][!(coefficient %like% 'ligation')] %>% pivot_wider(names_from = model_type, values_from = log_value) %>% as.data.table()

# temp
if ('Coefficients from model trained with EM\n(Unknown annotations)' %in% colnames(results_wide)) {
    results_wide[is.na(get('Coefficients from model trained with EM\n(Unknown annotations)')), paste('Coefficients from model trained with EM\n(Unknown annotations)') := Inf]
} else {
    results_wide[['Coefficients from model trained with EM\n(Unknown annotations)']] = Inf
}

colors <- brewer.pal(4, "Set2")
plot = ggplot(results_wide, aes(x = get('Coefficients from model trained without EM\n(Known annotations)'), y = get('Coefficients from model trained with EM\n(Unknown annotations)'), color = factor(simulator_type), shape = sequence_type)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray60", linewidth=1) +
  geom_point(size = 7) +
  scale_color_manual(values = colors) +
  labs(color = "MH simulation type", shape = "Sequence/Productivity type", x = 'Coefficients from model trained without EM\n(Known annotations)', y = 'Coefficients from model trained with EM\n(Unknown annotations)') +
  theme_cowplot(font_family = 'Arial') + 
  background_grid(major = 'xy') + 
  panel_border(color = 'gray60', size = 1.1) +
  theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 18), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), panel.spacing = unit(3, "lines"), legend.key.height = unit(2, 'cm'))  

# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_simulator_experiments/compare_EM_noEM_simulator/compare_trim_motif.pdf')

ggsave(file_name, plot = plot, width = 14, height = 8, units = 'in', dpi = 750, device = cairo_pdf)
