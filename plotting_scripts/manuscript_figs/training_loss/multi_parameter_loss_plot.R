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

PARAM_GROUP <<- 'nonproductive_v-j_trim'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))

NCPU <<- as.numeric(2)

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(1)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(2)

source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/model_evaluation_functions.R'))

model_types = c('motif_two-side-base-count-beyond_mh', 'motif_two-side-base-count-beyond', 'null')

# get all loss results
eval_results = data.table()
for (m in model_types) {
    if (m == 'null'){
        next
    } else if (m %like% 'mh'){
        l = 'True'
    } else {
        l = 'False'
    }
    path = get_model_eval_file_path(L2=l, model_type=m)
    temp = fread(path)
    temp$L2 = l
    eval_results = rbind(eval_results, temp)
}

# get null model loss results
LEFT_NUC_MOTIF_COUNT <<- 0
RIGHT_NUC_MOTIF_COUNT <<- 0

source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

path = get_model_eval_file_path(L2='False', model_type='null')
temp = fread(path)
temp$L2 = 'False'
eval_results = rbind(eval_results, temp)

# get model types, and make them neat
new_model_types = c('1x2motif + two-side base-count beyond + MH\n(38 params)', '1x2motif + two-side base-count-beyond\n(24 params)', 'null (0 params)')

# reassign model names
eval_results$model_type = mapvalues(eval_results$model_type, from = model_types, to = new_model_types) 

# make model names neater and set palette colors
colors = set_color_palette(c(new_model_types), with_params = TRUE)

# order loss types
loss_order = c('Log loss on training data', 'Expected log loss across training data')
eval_results$loss_type = factor(eval_results$loss_type, levels = loss_order)

# write plot data
cols = c('model_type', 'loss_type', 'log_loss')
loss_data = eval_results[, ..cols]
fwrite(loss_data, paste0(PROJECT_PATH, '/plotting_scripts/manuscript_figs/training_loss/loss.tsv'), sep = '\t')

# create plot
plot = plot_model_evaluation_loss_paracoord(eval_results, loss_bound = NULL, color_palette = colors, expand_var = 1.2)

plot = plot + ylab('Expected per-sequence log loss\n')

# save plot
file_name = paste0(PROJECT_PATH, '/plotting_scripts/manuscript_figs/training_loss/loss_compare.pdf')

ggsave(file_name, plot = plot, width = 32, height = 18, units = 'in', dpi = 750, device = cairo_pdf)
