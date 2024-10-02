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
library(ggpattern)
library(forcats)

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

file_root = paste0(PROJECT_PATH, '/plotting_scripts/manuscript_figs/validation_loss/trimMH_ligationMH_model/')

model_types = list('motif_two-side-base-count-beyond_average-mh_ligation-mh' = c(1,2), 'motif_two-side-base-count-beyond' = c(1, 2), 'motif_left-base-count-beyond_average-mh_ligation-mh' = c(1,2), 'motif_right-base-count-beyond_average-mh_ligation-mh' = c(1,2), 'null' = c(0,0), 'two-side-base-count-beyond_average-mh_ligation-mh' = c(0,0))
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
    for (annot in c('igor_alpha', 'validation_igor_alpha', 'validation_adaptive_gamma')){
        if (annot == 'igor_alpha'){
            # load predictions
            path = get_model_predictions_file_path(L2=L2, model_type = m, pred_type = 'training') 
        } else {
            path = get_validation_predictions_file_path(L2=L2, validation_annotation=annot, model_type = m)
        }
        temp = fread(path)
        if ('new_type' %in% colnames(temp)) {
            temp = temp[new_type != '0']
        }
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
loss_types = list('igor_alpha' = 'TRA training dataset', 'validation_igor_alpha' = 'TRA testing dataset', 'validation_adaptive_gamma' = 'TRG testing dataset')
eval_results$long_loss_type = mapvalues(eval_results$validation_data_type, from = names(loss_types), to = loss_types)
predictions$long_loss_type = mapvalues(predictions$validation_annotation_type, from = names(loss_types), to = loss_types)

# order loss types
eval_results$long_loss_type = factor(eval_results$long_loss_type, levels = loss_types)
predictions$long_loss_type = factor(predictions$long_loss_type, levels = loss_types)

eval_results = eval_results[validation_data_type %in% names(loss_types)]
predictions = predictions[validation_annotation_type %in% names(loss_types)]

eval_results[long_loss_type %like% 'training' & validation_productivity == 'nonproductive', training_subset := "sequences used for training model"]
predictions[long_loss_type %like% 'training' & validation_productivity == 'nonproductive', training_subset := "sequences used for training model"]

# calculate difference in loss
cols = c('validation_productivity', 'long_loss_type', 'log_loss', 'model_type', 'training_subset')
tog_diff = eval_results[, ..cols] %>% pivot_wider(names_from = 'model_type', values_from = 'log_loss') %>% as.data.table()

tog_diff[, 'full_to_no_MH' := get("motif_two-side-base-count-beyond_average-mh_ligation-mh")-get("motif_two-side-base-count-beyond")]
tog_diff[, 'full_to_no_motif' := get("motif_two-side-base-count-beyond_average-mh_ligation-mh")-get("two-side-base-count-beyond_average-mh_ligation-mh")]
tog_diff[, 'full_to_no_left_bc' := get("motif_two-side-base-count-beyond_average-mh_ligation-mh") - get("motif_right-base-count-beyond_average-mh_ligation-mh")]
tog_diff[, 'full_to_no_right_bc' := get("motif_two-side-base-count-beyond_average-mh_ligation-mh") - get("motif_left-base-count-beyond_average-mh_ligation-mh")]
tog_diff[, 'full_to_null' := get("motif_two-side-base-count-beyond_average-mh_ligation-mh")-get("null")]

cols2 = c('validation_productivity', 'long_loss_type', 'training_subset', 'full_to_no_MH', 'full_to_no_motif', 'full_to_no_left_bc', 'full_to_no_right_bc', 'full_to_null')
tog_diff_long = tog_diff[, ..cols2] %>% pivot_longer(starts_with('full_to'), names_to = 'diff_type', values_to = 'diff_log_loss') %>% as.data.table() 

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
mae_per = pred_subset[, mean(error), by = .(v_gene, j_gene, validation_productivity, long_loss_type, model_type, training_subset)]
mae = mae_per[, mean(V1), by = .(validation_productivity, long_loss_type, model_type, training_subset)]

## get difference in MAE
cols = c('validation_productivity', 'long_loss_type', 'V1', 'model_type', 'training_subset')
mae_diff = mae[, ..cols] %>% pivot_wider(names_from = 'model_type', values_from = 'V1') %>% as.data.table()

mae_diff[, 'full_to_no_MH' := get("motif_two-side-base-count-beyond_average-mh_ligation-mh")-get("motif_two-side-base-count-beyond")]
mae_diff[, 'full_to_no_motif' := get("motif_two-side-base-count-beyond_average-mh_ligation-mh")-get("two-side-base-count-beyond_average-mh_ligation-mh")]
mae_diff[, 'full_to_no_left_bc' := get("motif_two-side-base-count-beyond_average-mh_ligation-mh") - get("motif_right-base-count-beyond_average-mh_ligation-mh")]
mae_diff[, 'full_to_no_right_bc' := get("motif_two-side-base-count-beyond_average-mh_ligation-mh") - get("motif_left-base-count-beyond_average-mh_ligation-mh")]
mae_diff[, 'full_to_null' := get("motif_two-side-base-count-beyond_average-mh_ligation-mh")-get("null")]

cols2 = c('validation_productivity', 'long_loss_type', 'training_subset', 'full_to_no_MH', 'full_to_no_motif', 'full_to_no_left_bc', 'full_to_no_right_bc', 'full_to_null')

mae_diff_long = mae_diff[, ..cols2] %>% pivot_longer(starts_with('full_to'), names_to = 'diff_type', values_to = 'diff_mae') %>% as.data.table() 

# plot loss improvement
colormap = c('full_to_no_left_bc'='#66a61e', 'full_to_no_motif' = '#e7298a', 'full_to_no_right_bc'='#666666', 'full_to_no_MH'='#e6ab02')

tog_diff_long2 <- tog_diff_long[!(diff_type %like% 'null')] %>%
  group_by(validation_productivity, long_loss_type) %>%
  mutate(diff_type2 = factor(diff_type, levels = unique(diff_type))) %>%
  ungroup()

bar = ggplot(tog_diff_long2)+
      facet_grid(rows = vars(validation_productivity), cols = vars(long_loss_type), scales = "free_x", space = "free")+
      geom_col_pattern(aes(x = diff_type2, y = -diff_log_loss, pattern = training_subset, fill = diff_type, color = ifelse(training_subset == "sequences used for training model", "black", NA)),
                   size = 0.5,
                   pattern_density = 0.05,  # Adjust the density of the pattern
                   pattern_spacing = 0.05,  # Adjust the spacing of the pattern
                   pattern_fill = "black",  # Color of the pattern lines
                   pattern_color = "black",
                   pattern_size = 0.5,
                   pattern_angle = 45) +    # Angle of the stripes
      scale_pattern_manual(values = c("sequences used for training model" = "stripe", "NA" = "none"))+
      scale_color_identity()+
      scale_fill_manual(values = colormap, name = 'Model comparison') +
      theme_cowplot(font_family = 'Arial') + 
      xlab('') +
      ylab('Improvement in\nexpected per-sequence log loss') +
      background_grid(major = 'xy') + 
      panel_border(color = 'gray60', size = 1.5) +
      theme(text = element_text(size = 12), axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(size = 10), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 12), plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), legend.position = 'bottom', legend.direction = 'horizontal') 

ggsave(paste0(file_root, '/loss_dropout.pdf'), plot = bar, width = 8.4, height = 4, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)

mae_diff_long2 <- mae_diff_long[!(diff_type %like% 'null')] %>%
  group_by(validation_productivity, long_loss_type) %>%
  mutate(diff_type2 = factor(diff_type, levels = unique(diff_type))) %>%
  ungroup()

mae_bar = ggplot(mae_diff_long2)+
      facet_grid(rows = vars(validation_productivity), cols = vars(long_loss_type), scales = "free_x", space = "free")+
      geom_col_pattern(aes(x = diff_type2, y = -diff_mae, pattern = training_subset, fill = diff_type, color = ifelse(training_subset == "sequences used for training model", "black", NA)),
                   size = 0.5,
                   pattern_density = 0.05,  # Adjust the density of the pattern
                   pattern_spacing = 0.05,  # Adjust the spacing of the pattern
                   pattern_fill = "black",  # Color of the pattern lines
                   pattern_color = "black",
                   pattern_size = 0.5,
                   pattern_angle = 45) +    # Angle of the stripes
      scale_pattern_manual(values = c("sequences used for training model" = "stripe", "NA" = "none"))+
      scale_color_identity()+
      scale_fill_manual(values = colormap, name = 'Model comparison') +
      theme_cowplot(font_family = 'Arial') + 
      xlab('') +
      ylab('Improvement in MAE')+
      background_grid(major = 'xy') + 
      panel_border(color = 'gray60', size = 1.5) +
      theme(text = element_text(size = 12), axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(size = 10), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 12), plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), legend.position = 'bottom', legend.direction = 'horizontal') 

ggsave(paste0(file_root, '/mae_dropout.pdf'), plot = mae_bar, width = 8.4, height = 4, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)

