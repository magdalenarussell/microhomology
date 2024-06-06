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
TRIMMING_PROB_MODEL <<- c('uniform')
SAMPLE = c(FALSE)

results = data.table()

for (param_r in PARAM_root){
    PARAM_GROUP <<- paste0(param_r, '_v-j_trim_ligation-mh')
    source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))
    for (s in SAMPLE){
        for (trim_model in TRIMMING_PROB_MODEL){
            for (param in LIGATION_PARAMS){
                for (param2 in INT_MH_PARAMS){

                    ANNOTATION_TYPE <<- 'igor_mh_sim_alpha'
                    old_annotation = ANNOTATION_TYPE

                    source(paste0(MOD_PROJECT_PATH, '/config/file_paths.R'))
                    source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

                    ANNOTATION_TYPE <<- paste0(old_annotation, '_from_', trim_model, '_MHprob', param, '_trimMHprob', param2)
                    # Read in model coefficient data 
                    coef_path = get_model_coef_file_path(L2, sample_annotation=s)
                    if (!file.exists(coef_path)){
                        assign(paste0('together', param), NULL)

                        next
                    }
                    coefs = fread(coef_path)
                    temp = coefs[coefficient %like% 'mh']
                    temp$sequence_type = param_r
                    temp$annotation_sampled = s
                    temp$trimming_model = trim_model
                    temp$lig_MH_param = param
                    temp$trim_MH_param = param2

                    results = rbind(results, temp)
                }
            }
        }
    }
}

results$log_value = results$value/log(10)
cols = c('coefficient', 'log_value', 'sequence_type', 'trimming_model', 'lig_MH_param', 'trim_MH_param', 'annotation_sampled')
results_wide = results[, ..cols][coefficient %like% 'ligation'] %>% pivot_wider(names_from = sequence_type, values_from = log_value) %>% as.data.table()

results_wide[, MH_params := paste0('ligMH param = ', lig_MH_param, ', trimMH param = ', trim_MH_param)]

colors <- brewer.pal(4, "Set2")
plot = ggplot(results_wide, aes(x = both, y = nonproductive, color = factor(MH_params))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray60", linewidth=1) +
  geom_point(size = 7) +
  scale_color_manual(values = colors) +
  labs(x = "ligMH parameters from model trained with\nall sequences", y = "ligMH parameters from model trained with\nnonproductive sequences", color = "MH param", shape = "Sequence type") +
  theme_cowplot(font_family = 'Arial') + 
  background_grid(major = 'xy') + 
  panel_border(color = 'gray60', size = 1.1) +
  theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 18), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), panel.spacing = unit(3, "lines"), legend.key.height = unit(2, 'cm'))  

# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_simulator_experiments/motif_base_count_mh-config-count_ligation-mh_trims_from_uniform_experiment_EM_ligMH_trimMH_sims/compare_lig.pdf')

ggsave(file_name, plot = plot, width = 14, height = 8, units = 'in', dpi = 750, device = cairo_pdf)

cols = c(cols, 'position')
results_wide = results[, ..cols][!(coefficient %like% 'ligation')] %>% pivot_wider(names_from = sequence_type, values_from = log_value) %>% as.data.table()

results_wide[, MH_params := paste0('ligMH param = ', lig_MH_param, ', trimMH param = ', trim_MH_param)]

colors <- brewer.pal(4, "Set2")
plot = ggplot(results_wide, aes(x = both, y = nonproductive, color = factor(MH_params), shape = position)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray60", linewidth=1) +
  geom_point(size = 7) +
  scale_color_manual(values = colors) +
  labs(x = "intMH parameters from model trained with\nall sequences", y = "intMH parameters from model trained with\nnonproductive sequences", color = "MH param", shape = "Overlap region") +
  theme_cowplot(font_family = 'Arial') + 
  background_grid(major = 'xy') + 
  panel_border(color = 'gray60', size = 1.1) +
  theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 18), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), panel.spacing = unit(3, "lines"), legend.key.height = unit(2, 'cm'))  

# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_simulator_experiments/motif_base_count_mh-config-count_ligation-mh_trims_from_uniform_experiment_EM_ligMH_trimMH_sims/compare_int.pdf')

ggsave(file_name, plot = plot, width = 14, height = 8, units = 'in', dpi = 750, device = cairo_pdf)


