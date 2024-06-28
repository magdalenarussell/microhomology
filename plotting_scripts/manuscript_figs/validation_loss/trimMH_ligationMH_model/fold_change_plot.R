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

ANNOTATION_TYPE <<- 'igor_alpha' 

PARAM_GROUP <<- 'nonproductive_v-j_trim_ligation-mh'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))

NCPU <<- as.numeric(2)

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(1)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(2)
MODEL_TYPE = 'motif_two-side-base-count-beyond_average-mh_ligation-mh' 

source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

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
eval_results[model_type == '', model_type := 'null']

# get predictions
predictions = data.table()
for (m in names(model_types)){
    LEFT_NUC_MOTIF_COUNT <<- model_types[[m]][1]
    RIGHT_NUC_MOTIF_COUNT <<- model_types[[m]][2]
    MODEL_TYPE <<- m
    for (annot in c('igor_alpha', 'validation_adaptive_alpha', 'validation_igor_alpha')){
        if (annot == 'igor_alpha'){
            # load predictions
            path = get_model_predictions_file_path(L2=L2, model_type = m) 
        } else {
            path = get_validation_predictions_file_path(L2=L2, validation_annotation=annot, model_type = m)
        }
        temp = fread(path)[new_type != '0']
        cols = c('v_gene', 'j_gene', 'v_trim', 'j_trim', 'ligation_mh', 'index', 'predicted_prob')
        temp = temp[, ..cols]
        for (prod in c('nonproductive_v-j_trim_ligation-mh', 'productive_v-j_trim_ligation-mh')) {
            # load empirical data
            old_annot <<- 'igor_alpha'
            old_prod <<- 'nonproductive_v-j_trim_ligation-mh'
            ANNOTATION_TYPE <<- annot
            PARAM_GROUP <<- prod
            path = processed_data_path()
            temp2 = fread(path)
            tog = merge(temp2, temp)
            tog$validation_productivity = str_split(prod, '_')[[1]][1] 
            tog$validation_param_group = prod
            tog$validation_annotation_type = annot
            tog$model_type = m
            predictions = rbind(predictions, tog)
            ANNOTATION_TYPE <<- old_annot
            PARAM_GROUP <<- old_prod
        }
    }
}

# reassign loss types
loss_types = list('igor_alpha' = 'TRA training\ndataset\n(IGoR)', 'validation_igor_alpha' = 'TRA testing\ndataset\n(IGoR)', 'validation_adaptive_alpha' = 'TRA testing\ndataset\n(Adaptive)')
eval_results$long_loss_type = mapvalues(eval_results$validation_data_type, from = names(loss_types), to = loss_types)
predictions$long_loss_type = mapvalues(predictions$validation_annotation_type, from = names(loss_types), to = loss_types)

# order loss types
eval_results$long_loss_type = factor(eval_results$long_loss_type, levels = loss_types)
predictions$long_loss_type = factor(predictions$long_loss_type, levels = loss_types)

# calculate fold change in loss
cols = c('validation_productivity', 'long_loss_type', 'log_loss', 'model_type')
tog_fold = eval_results[, ..cols] %>% pivot_wider(names_from = 'model_type', values_from = 'log_loss') %>% as.data.table()
tog_fold[, paste0('relative to motif + base-count model') := get("motif_two-side-base-count-beyond_average-mh_ligation-mh")/get("motif_two-side-base-count-beyond")]
tog_fold[, paste0('relative to null model') := get("motif_two-side-base-count-beyond_average-mh_ligation-mh")/get("null")]

cols2 = c('validation_productivity', 'long_loss_type', 'relative to motif + base-count model', 'relative to null model')
tog_fold_long = tog_fold[, ..cols2] %>% pivot_longer(starts_with('relative'), names_to = 'fold_type', values_to = 'fold_loss_change') %>% as.data.table() 

# calculate MAE
pred_subset = predictions[!(is.na(count)) & !(is.na(weighted_observation))]
## normalize predicted probabilities
pred_subset[, pred_sum := sum(predicted_prob), by = .(v_gene, j_gene, validation_productivity, long_loss_type, model_type)]
pred_subset[, norm_predicted_prob := predicted_prob/pred_sum]
## get index weights
pred_subset[, index_sum := sum(norm_predicted_prob), by = .(index, validation_productivity, long_loss_type, model_type)]
pred_subset[, index_weight := norm_predicted_prob/index_sum]
## normalize empirical probabilities
pred_subset[, weighted_count := count * index_weight]
pred_subset[, emp_sum := sum(weighted_count), by = .(v_gene, j_gene, validation_productivity, long_loss_type, model_type)]
pred_subset[, norm_empirical_prob := weighted_count/emp_sum]

## get per-gene-pair MAE
pred_subset[, error := abs(norm_predicted_prob - norm_empirical_prob)]
mae_per = pred_subset[, mean(error), by = .(v_gene, j_gene, validation_productivity, long_loss_type, model_type)]
mae = mae_per[, mean(V1), by = .(validation_productivity, long_loss_type, model_type)]

## get fold MAE
cols = c('validation_productivity', 'long_loss_type', 'V1', 'model_type')
mae_fold = mae[, ..cols] %>% pivot_wider(names_from = 'model_type', values_from = 'V1') %>% as.data.table()
mae_fold[, paste0('relative to motif + base-count model') := get("motif_two-side-base-count-beyond_average-mh_ligation-mh")-get("motif_two-side-base-count-beyond")]
mae_fold[, paste0('relative to null model') := get("motif_two-side-base-count-beyond_average-mh_ligation-mh")-get("null")]

cols2 = c('validation_productivity', 'long_loss_type', 'relative to motif + base-count model', 'relative to null model')
mae_fold_long = mae_fold[, ..cols2] %>% pivot_longer(starts_with('relative'), names_to = 'fold_type', values_to = 'fold_mae_change') %>% as.data.table() 

for (prod in c('nonproductive', 'productive')){
    plot_loss = ggplot(tog_fold_long[validation_productivity == prod]) +
             geom_point(aes(x = long_loss_type, y = fold_loss_change, color = fold_type), size = 6) +
             geom_line(aes(x = long_loss_type, y = fold_loss_change, color = fold_type, group = fold_type), linewidth = 3) +
             geom_hline(yintercept = 1, linetype = 'dashed', linewidth = 1, color = 'black') +
             theme_cowplot(font_family = 'Arial') + 
             xlab(' ') +
             ylab('Fold change in\nexpected per-sequence log loss\n') +
             background_grid(major = 'xy') + 
             panel_border(color = 'gray60', size = 1.5) +
             theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 18), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
             scale_color_brewer(palette = 'Dark2') +
             labs(color = 'Fold change type')

    file_name = paste0(PROJECT_PATH, '/plotting_scripts/manuscript_figs/validation_loss/trimMH_ligationMH_model/loss_fold_compare_', prod, '.pdf')

    ggsave(file_name, plot = plot_loss, width = 14, height = 7, units = 'in', dpi = 750, device = cairo_pdf)
}

for (prod in c('nonproductive', 'productive')){
    plot_mae = ggplot(mae_fold_long[validation_productivity == prod]) +
             geom_point(aes(x = long_loss_type, y = fold_mae_change, color = fold_type), size = 6) +
             geom_line(aes(x = long_loss_type, y = fold_mae_change, color = fold_type, group = fold_type), linewidth = 3) +
             geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 1, color = 'black') +
             theme_cowplot(font_family = 'Arial') + 
             xlab(' ') +
             ylab('Difference in MAE\n') +
             background_grid(major = 'xy') + 
             panel_border(color = 'gray60', size = 1.5) +
             theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 18), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
             scale_color_brewer(palette = 'Dark2') +
             labs(color = 'Difference type')

    file_name = paste0(PROJECT_PATH, '/plotting_scripts/manuscript_figs/validation_loss/trimMH_ligationMH_model/mae_fold_compare_', prod, '.pdf')

    ggsave(file_name, plot = plot_mae, width = 14, height = 7, units = 'in', dpi = 750, device = cairo_pdf)
}
