source('config/config.R')

library(cli)
# library(devtools)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(foreach)
library(doParallel)
library(tidyverse)
library(plyr)
library(data.table)
setDTthreads(1)
library(Biostrings)
# library(mclogit)
# library(matrixcalc)
# library(RhpcBLASctl)
# omp_set_num_threads(1)
# blas_set_num_threads(1)

ANNOTATION_TYPE <<- 'igor_alpha' 

PARAM_GROUP <<- 'nonproductive_v-j_trim_ligation-mh'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))

NCPU <<- as.numeric(2)

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(1)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(2)

source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

model_types = list('motif_two-side-base-count-beyond_average-mh_ligation-mh' = c(1,2), 'motif_two-side-base-count-beyond' = c(1, 2), 'null' = c(0,0))

L2 = 'True'

# get all loss results
eval_results = data.table()
for (m in names(model_types)) {
    LEFT_NUC_MOTIF_COUNT <<- model_types[[m]][1]
    RIGHT_NUC_MOTIF_COUNT <<- model_types[[m]][2]
    path = get_model_eval_file_path(L2=L2, model_type=m)
    temp = fread(path)
    temp$L2 = L2
    eval_results = rbind(eval_results, temp)
}

# get model types, and make them neat
# reassign model names
new_model_types = list('motif_two-side-base-count-beyond_average-mh_ligation-mh' = '1x2motif + two-side base-count beyond +\naverage trim MH + ligation MH\n(26 params)', 'motif_two-side-base-count-beyond' = '1x2motif + two-side base-count-beyond\n(24 params)', 'null' = 'null (0 params)')
eval_results$long_model_type = mapvalues(eval_results$model_type, from = names(new_model_types), to = unlist(new_model_types))
eval_results[model_type == '', long_model_type := 'null (0 params)']

# reassign loss types
loss_types = list('igor_alpha' = 'TRA training dataset\n(IGoR-annotated)', 'validation_adaptive_alpha' = 'TRA testing dataset\n(Adaptive-annotated)')
eval_results$long_loss_type = mapvalues(eval_results$validation_data_type, from = names(loss_types), to = loss_types)

# order loss types
eval_results$long_loss_type = factor(eval_results$long_loss_type, levels = loss_types)

# write plot data
cols = c('long_model_type', 'long_loss_type', 'log_loss')
loss_data = eval_results[, ..cols]
fwrite(loss_data, paste0(PROJECT_PATH, '/plotting_scripts/manuscript_figs/validation_loss/trimMH_ligationMH_model/loss.tsv'), sep = '\t')

# create plot
colors = set_color_palette(unlist(new_model_types), with_params = TRUE)
plot_NP = plot_model_evaluation_loss_paracoord(eval_results[validation_productivity == 'nonproductive'], loss_bound = NULL, expand_var = 1.2, color_palette = colors, productivity_facet = FALSE)

plot_NP = plot_NP + ylab('Expected per-sequence log loss\n')

# add reference line
subset = eval_results[validation_productivity == 'nonproductive'][long_loss_type %like% 'train'][long_model_type %like% 'MH']

plot_NP = plot_NP + geom_hline(data = subset, aes(yintercept = log_loss, color = long_model_type), linewidth = 2, linetype = 'dashed')

# save plot
file_name = paste0(PROJECT_PATH, '/plotting_scripts/manuscript_figs/validation_loss/trimMH_ligationMH_model/loss_compare_nonproductive.pdf')

ggsave(file_name, plot = plot_NP, width = 32, height = 18, units = 'in', dpi = 750, device = cairo_pdf)

plot_P = plot_model_evaluation_loss_paracoord(eval_results[validation_productivity == 'productive'], loss_bound = NULL, expand_var = 1.2, color_palette = colors, productivity_facet = FALSE)

plot_P = plot_P + ylab('Expected per-sequence log loss\n')

plot_P = plot_P + geom_hline(data = subset, aes(yintercept = log_loss, color = long_model_type), linewidth = 2, linetype = 'dashed')

# save plot
file_name = paste0(PROJECT_PATH, '/plotting_scripts/manuscript_figs/validation_loss/trimMH_ligationMH_model/loss_compare_productive.pdf')

ggsave(file_name, plot = plot_P, width = 32, height = 18, units = 'in', dpi = 750, device = cairo_pdf)
