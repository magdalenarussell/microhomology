source('config/config.R')

library(foreach)
library(doParallel)
library(tidyverse)
library(plyr)
library(data.table)
setDTthreads(1)
library(Biostrings)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(mclogit)
library(matrixcalc)
library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)

ANNOTATION_TYPE <<- 'igor' 

TRIM_TYPE <<- 'v_trim'

PRODUCTIVITY <<- 'nonproductive' 

MOTIF_TYPE <<- 'unbounded' 

NCPU <<- 2

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
stopifnot(GENE_NAME == 'v_gene')

MODEL_GROUP <<- 'all_subjects'

GENE_WEIGHT_TYPE <<- 'p_gene_marginal' 

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 0
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 0 

UPPER_TRIM_BOUND <<- as.numeric(14) 

MODEL_TYPE <<- 'two-side-base-count'

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(10)

TYPE <<- 'v_gene_family_loss'
source(paste0(MOD_PROJECT_PATH, 'scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'scripts/model_fitting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'plotting_scripts/plotting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'plotting_scripts/residual_comparison_functions.R'))
source(paste0(MOD_PROJECT_PATH,'scripts/analysis_scripts/pwm_profile_functions.R'))
source(paste0(MOD_PROJECT_PATH,'scripts/model_evaluation_functions.R'))

genes = c('TRBV9', 'TRBV13')

# Read in dist and residual data for the two-side-base-count model
predicted_trims1 = get_predicted_distribution_data() 
per_gene_resid1 = calculate_rmse_by_gene(predicted_trims1) 
per_gene_resid1$model = MODEL_TYPE
per_gene_trim_resid1 = calculate_rmse_by_gene_trim(predicted_trims1)
per_gene_trim_resid1$model = MODEL_TYPE
predicted_trims1$model = 'two-side-base-count model,\t\nall training data' 

# set parameters and re-source files for the motif + two-side-base-count-beyond model
MODEL_TYPE1 = MODEL_TYPE
LEFT_SIDE1 = LEFT_SIDE_TERMINAL_MELT_LENGTH
MODEL_TYPE <<- 'motif_two-side-base-count-beyond'
stopifnot(MODEL_TYPE %like% 'motif')
# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(1)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(2)

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('plotting_scripts/plotting_functions.R')
source('scripts/analysis_scripts/pwm_profile_functions.R')

# Read in dist and residual data for the motif + two-side-base-count-beyond model
predicted_trims2 = get_predicted_distribution_data() 
per_gene_resid2 = calculate_rmse_by_gene(predicted_trims2) 
per_gene_resid2$model = MODEL_TYPE
per_gene_trim_resid2 = calculate_rmse_by_gene_trim(predicted_trims2)
per_gene_trim_resid2$model = MODEL_TYPE
predicted_trims2$model = '+ 1x2 motif,\nall training data\t' 

#merge two predictions for both models
together = rbind(per_gene_resid1, per_gene_resid2)
together_trim = rbind(per_gene_trim_resid1, per_gene_trim_resid2)
predictions = rbind(predicted_trims1, predicted_trims2, fill = TRUE)

# order observations based on model type
together$model = factor(together$model, levels = c(MODEL_TYPE1, MODEL_TYPE))

# convert data to wide format for plotting
wide = together[, -c('p_gene')] %>% pivot_wider(names_from = 'model', values_from = 'rmse') %>% as.data.table()

# determine outlier threshold
slopes = wide[,(get(MODEL_TYPE) - get(MODEL_TYPE1)), by = gene]
setnames(slopes, 'V1', 'slope')
together = merge(together, slopes)
together = together[order(abs(slope))]
outlier_count = ceiling(nrow(together)*0.1)
if ((outlier_count %% 2) != 0){
    outlier_count = outlier_count + 1
}
outliers = together[(nrow(together)-outlier_count + 1):nrow(together)]
no_outliers = together[1:(nrow(together)-outlier_count)]
cutoff = mean(c(min(no_outliers$slope), max(outliers$slope)))

# calculate the difference
wide$diff = wide[['motif_two-side-base-count-beyond']] - wide[['two-side-base-count']]

fwrite(wide, paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/motif_exploration/hist.tsv'), sep = '\t')

# create histogram plot
plot = ggplot(wide) +
    geom_histogram(aes(x = diff), binwidth= 0.01) +
    geom_vline(xintercept = cutoff, size = 5, color = 'gray') +
    geom_text(x = -0.15, y = 6, label = '\"improved\ngenes\"', color = 'gray', size = 12, lineheight = 0.8, check_overlap = TRUE)+
    ylab('Gene count\n') +
    xlab('\nPer-gene RMSE difference') +
    theme_cowplot(font_family = 'Arial') +
    theme(text = element_text(size = 30), axis.text.x=element_text(size = 30), axis.text.y = element_text(size = 30), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + 
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) 

# re-fit_model without outliers
motif_data = aggregate_all_subject_data()
new_motif_data = motif_data[gene %in% no_outliers$gene]

# plot predicted distributions on held out (outlier genes)
model = fit_model(new_motif_data)
outlier_motif_data = motif_data[gene %in% outliers$gene]
outlier_motif_data = process_data_for_model_fit(outlier_motif_data)
outlier_motif_data$predicted_prob = temp_predict(model, newdata = outlier_motif_data)
outlier_motif_data[, empirical_prob := count/sum(count), by = .(subject, gene)]
outlier_motif_data$model = paste0('+ 1x2motif,\nremoving \"improved genes\"\t')

predictions = rbind(predictions, outlier_motif_data, fill = TRUE)

# re-fit_model without outliers and genes similar to outliers
gene_families = get_gene_families(cluster_count = 5, combine_by_terminal = FALSE, full_sequence = TRUE, align = TRUE)$cluster_data
outlier_clusters = unique(gene_families[gene %in% unique(outliers$gene)]$clusters_grouped)
outlier_cluster_genes = unique(gene_families[clusters_grouped %in% outlier_clusters]$gene)
no_outliers = together[!(gene %in% outliers) & !(gene %in% outlier_cluster_genes)]
new_motif_data = motif_data[gene %in% no_outliers$gene]

# plot predicted distributions on held out (outlier genes + similar)
model = fit_model(new_motif_data)
outlier_motif_data = motif_data[gene %in% outliers$gene]
outlier_motif_data = process_data_for_model_fit(outlier_motif_data)
outlier_motif_data$predicted_prob = temp_predict(model, newdata = outlier_motif_data)
outlier_motif_data[, empirical_prob := count/sum(count), by = .(subject, gene)]
outlier_motif_data$model = paste0('+ 1x2motif,\nremoving \"improved genes\" + similar')

predictions = rbind(predictions, outlier_motif_data, fill = TRUE)

# order model types
model_types = c('two-side-base-count model,\t\nall training data', '+ 1x2motif,\nremoving \"improved genes\" + similar', '+ 1x2motif,\nremoving \"improved genes\"\t', '+ 1x2 motif,\nall training data\t')
predictions$model = factor(predictions$model, levels = model_types) 

# set colors for each model
colors = c('#8c510a','#80cdc1', '#35978f', '#01665e')
names(colors) = model_types

cols = c('model', 'subject', 'gene', 'trim_length', 'empirical_prob', 'predicted_prob')
plot_data = unique(predictions[gene %in% genes][, ..cols])
fwrite(plot_data, paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/motif_exploration/dists.tsv'), sep = '\t')

# plot distributions for each gene
for (gene_name in genes){
    # extract only gene-specific data
    important_cols = c('trim_length', 'predicted_prob', 'gene', 'model')
    predicted_data = unique(predictions[gene == gene_name, ..important_cols])
    predicted_data$model = factor(predicted_data$model, levels = model_types)
    empirical_data = predictions[gene == gene_name][order(subject, trim_length)]
    empirical_data$model = factor(empirical_data$model, levels = model_types)
    title = paste0(gene_name)

    # get gene sequence
    gene_seq = get_gene_sequence(gene_name, max(predictions$trim_length))
    gene_seq_with_positions = get_plot_positions_for_gene_sequence(gene_seq)
    gene_seq_with_positions =get_motif_colors(gene_seq_with_positions, c('CTT', 'CGT'), 'motif')
   
    max_prob = max(max(empirical_data$empirical_prob), max(predicted_data$predicted_prob))

    labels = data.table(model = model_types, yvar = max_prob, leftx = -2.1, rightx = UPPER_TRIM_BOUND + 0.1) 

    temp_plot = ggplot() +
        geom_line(data = empirical_data, aes(x = trim_length, y = empirical_prob, group = subject), size = 1.6, alpha = 0.3, color = 'gray60') +
        geom_line(data = predicted_data, aes(x = trim_length, y = predicted_prob, group = model, color = model), size = 2.75, alpha = 0.9) +
        # facet_grid(cols = vars(factor(model, levels = model_types))) +
        geom_vline(xintercept = 0, color = 'black', size = 3) +
        geom_text(data = gene_seq_with_positions, y = max_prob, aes(x = position, label = base, color = color), size = 10) +
        geom_text(data = labels, aes(y = yvar, x = leftx), label = '3\'- ', size = 6) +
        geom_text(data = labels, aes(y = yvar, x = rightx), label = ' -5\'', size = 6) +
        ggtitle(title) +
        xlab('Number of trimmed nucleotides') +
        ylab('Probability\n') +
        theme_cowplot(font_family = 'Arial') + 
        theme(text = element_text(size = 30), axis.text.x=element_text(size = 25), axis.text.y = element_text(size = 25), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5), legend.direction = 'horizontal', legend.position = 'bottom') + 
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) +
        ylim(c(0, max_prob + 0.02))+
        scale_color_manual(values = c(colors, motif = '#E79737', black = 'black'), breaks = c('two-side-base-count model,\t\nall training data', '+ 1x2 motif,\nall training data\t', '+ 1x2motif,\nremoving \"improved genes\"\t', '+ 1x2motif,\nremoving \"improved genes\" + similar')) +
        labs(color = 'Model,\ntraining data set\t')+
        guides(color = guide_legend(override.aes = list(size = 8))) 
    assign(gene_name, temp_plot)
} 

# isolate legend
legend = get_legend(TRBV9)
TRBV9 = TRBV9 + theme(legend.position = 'none')
TRBV13 = TRBV13 + theme(legend.position = 'none')

# align and combine all plots in grid
all = align_plots(TRBV9, TRBV13, plot, align = 'h', axis = 'lrtb')
first_row = plot_grid(all[[3]], nrow = 1, labels = c('A', ''), label_size = 40)
second_row = plot_grid(all[[1]], NULL, all[[2]], nrow = 1, rel_widths = c(1, 0.02, 1), labels = c('B', '', 'C'), label_size = 40)
second_row2 = plot_grid(second_row, NULL, legend, NULL, nrow = 4, rel_heights = c(1, 0.02, 0.1, 0.02), align = 'v', axis = 'rl')
all_tog = plot_grid(first_row, NULL, second_row2, nrow = 3, rel_heights = c(0.6, 0.02, 0.75))

# save grid plot
path = get_manuscript_path()
file_name = paste0(path, '/motif_exploration.pdf')
ggsave(file_name, plot = all_tog, width = 27, height = 15, units = 'in', dpi = 750, device = cairo_pdf)
