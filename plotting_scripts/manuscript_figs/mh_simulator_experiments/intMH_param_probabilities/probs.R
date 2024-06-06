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

PARAM_GROUP <<- 'nonproductive_v-j_trim'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))

NCPU <<- as.numeric(2)

L2 <<- 'True'

INT_MH_PARAMS <<- c(0, 0.1, 1, 10)


source(paste0(MOD_PROJECT_PATH,'/mh_simulation_scripts/intermediate-mh_signal_simulator_scripts/intermediate-mh_simulator_functions.R'))

probs = data.table()

for (mh in INT_MH_PARAMS){
    temp_df = data.table(intermediate_mh = seq(0, 10), covariate_func = mh * seq(0, 10))
    temp_df$int_mh_param = mh
    probs = rbind(probs, temp_df)
}

probs[, rel_prob := exp(covariate_func)/sum(exp(covariate_func)), by = int_mh_param]

probs[, zero_rel_prob := min(rel_prob), by = int_mh_param]
probs[, rel_prob_effect := rel_prob/zero_rel_prob]
probs[int_mh_param == 0, int_mh_param := 0.01]

# scatter = ggplot(probs)+
#     geom_point(aes(x = intermediate_mh, y = log10(rel_prob_effect), color = int_mh_param), size = 10) +
#     geom_line(aes(x = intermediate_mh, y = log10(rel_prob_effect), group = int_mh_param, color = int_mh_param), size = 5) +
#     xlab('Amount of intermediate MH') +
#     ylab('log10(relative change in trimming probability)') +
#     theme_cowplot(font_family = 'Arial') + 
#     labs(color = 'intMH simulation parameter') +
#     background_grid(major = 'xy') + 
#     scale_color_distiller(palette = 'Oranges', direction = +1, trans = "log10", oob = scales::oob_squish, limits = c(1e-2, 10)) +
#     panel_border(color = 'gray60', size = 1.1) +
#     theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 16), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), panel.spacing = unit(3, "lines"), legend.key.height = unit(2, 'cm'), legend.key.width = unit(1.5, 'cm'))  

# # save plot
# file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_simulator_experiments/intMH_param_probabilities/params_log.pdf')

ggsave(file_name, plot = scatter, width = 15, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

scatter = ggplot(probs)+
    geom_point(aes(x = intermediate_mh, y = covariate_func, color = int_mh_param), size = 10) +
    geom_line(aes(x = intermediate_mh, y = covariate_func, group = int_mh_param, color = int_mh_param), size = 5) +
    xlab('Amount of intermediate MH') +
    ylab('Relative covariate function value for\nintermediate-MH trimming probability adjustment') +
    theme_cowplot(font_family = 'Arial') + 
    labs(color = 'intMH simulation parameter') +
    background_grid(major = 'xy') + 
    scale_color_distiller(palette = 'Oranges', direction = +1, trans = "log10", oob = scales::oob_squish, limits = c(1e-2, 10)) +
    panel_border(color = 'gray60', size = 1.1) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 16), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), panel.spacing = unit(3, "lines"), legend.key.height = unit(2, 'cm'), legend.key.width = unit(1.5, 'cm'))  

# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_simulator_experiments/intMH_param_probabilities/params.pdf')

ggsave(file_name, plot = scatter, width = 15, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

