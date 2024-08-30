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

PARAM_GROUP <<- 'nonproductive_v-j_trim_ligation-mh'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))

NCPU <<- as.numeric(2)

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(1)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(2)

MODEL_TYPE <<- 'motif_two-side-base-count-beyond_average-mh_ligation-mh'

L2 <<- 'True'

ANNOTATION_TYPE <<- 'igor_alpha'

source(paste0(MOD_PROJECT_PATH, '/config/file_paths.R'))
source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

file_root = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/motif_base_count_average-mh_ligation-mh_model/', PARAM_GROUP)

# Read in model coefficient data 
coef_path = get_model_coef_file_path(L2)
coefs = fread(coef_path)
fwrite(coefs, paste0(file_root, '/coefs.tsv'), sep = '\t')

bounds = c(-0.450, 0.450)

map_positions_to_values_orient <- function(pwm, orient){
    if (orient == '5to3'){
        sign5 = -1
        sign3 = 1
    } else if (orient == '3to5'){
        sign5=1
        sign3=-1
    }
    pwm[side == '5end', position_value := sign5 * as.numeric(substring(position, nchar(position), nchar(position)))]
    pwm[side == '3end', position_value := sign3 *as.numeric(substring(position, nchar(position), nchar(position)))]
    return(pwm)
}

plot_oriented_motif_coefficient_heatmap <- function(model_coef_matrix, orient = '5to3', limits){
    model_coef_matrix = model_coef_matrix[coefficient == 'motif']
    model_coef_matrix = map_positions_to_values_orient(model_coef_matrix, orient)

    # convert to log_10
    model_coef_matrix$log_10_pdel = model_coef_matrix$value/log(10)
    # order variables
    model_coef_matrix$base = factor(model_coef_matrix$base, levels = c('T', 'G', 'C', 'A'))

    motif_length = RIGHT_NUC_MOTIF_COUNT + LEFT_NUC_MOTIF_COUNT

    if (orient == '5to3'){
        left_annot = "5\'"
        right_annot = "3\'"
        vline = LEFT_NUC_MOTIF_COUNT + 0.5 
        model_coef_matrix$position_value = factor(model_coef_matrix$position_value)
    } else if (orient == '3to5'){
        vline = RIGHT_NUC_MOTIF_COUNT + 0.5 
        left_annot = "3\'"
        right_annot = "5\'"
        posvals = unique(model_coef_matrix$position_value) 
        model_coef_matrix$position_value = factor(-1*model_coef_matrix$position_value, levels = -1*posvals[order(-(-1*posvals))])
    }
    
    plot = ggplot(model_coef_matrix, aes(x=position_value, y=base, fill=log_10_pdel)) +
        geom_tile() +
        theme_cowplot(font_family = 'Arial') + 
        xlab('Position') +
        ylab ('Base') +
        theme(text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 18), axis.title = element_text(size = 24))+
        geom_vline(xintercept = vline, size = 3.5, color = 'black') +
        guides(fill = guide_colourbar(barheight = 14)) +
        scale_fill_distiller(palette = 'PuOr', name = 'log10(probability of deletion)', limits = limits) +
        annotate("text", x = 0.35, y = 0.25, label = left_annot, size = 8) +  
        annotate("text", x = motif_length + 0.65, y = 0.25, label = right_annot, size = 8) +  
        coord_cartesian(ylim = c(1, 4), clip = "off")
    
    return(plot)
}

plot_indiv_base_count_coefficient <- function(model_coef_matrix, limits){
    model_coef_matrix = model_coef_matrix[(coefficient %like% 'base_count')]

    # convert to log_10
    model_coef_matrix$log_10_pdel = model_coef_matrix$value/log(10)
    
    # order variables
    unique_bases = c('AT', 'GC')
    left_bases = unique(model_coef_matrix[side %like% '5end']$base)
    right_bases = unique(model_coef_matrix[side %like% '3end']$base)
    trim_type = unique(model_coef_matrix$trim_type)

    if (length(left_bases) < 2){
        missing = unique_bases[!(unique_bases == left_bases)]
        fake_row = data.table(value = rep(NA, length(missing)), coefficient = rep("base_count", length(missing)), base = missing, position = rep('', length(missing)), side = rep('5end', length(missing)), log_10_pdel = rep(NA,length(missing)), trim_type = trim_type)
        model_coef_matrix = rbind(model_coef_matrix, fake_row, fill = TRUE)
    } 
    if (length(right_bases) < 2){
        missing = unique_bases[!(unique_bases == right_bases)]
        fake_row = data.table(value = rep(NA, length(missing)), coefficient = rep("base_count", length(missing)), base = missing, position = rep('', length(missing)), side = rep('3end', length(missing)), log_10_pdel = rep(NA,length(missing)), trim_type = trim_type)
        model_coef_matrix = rbind(model_coef_matrix, fake_row, fill = TRUE)
    } 
    
    model_coef_matrix[, long_name := paste0(substring(side, 1, 1), '\'-end base count')]

    for (side_long in unique(model_coef_matrix$long_name)){
        if (side_long %like% '5'){
            name = '5endplot'
        } else if (side_long %like% '3'){
            name = '3endplot'
        }
    
        plot = ggplot(model_coef_matrix[long_name == side_long], aes(x=long_name, y=base, fill=log_10_pdel)) +
            geom_tile() +
            theme_cowplot(font_family = 'Arial') + 
            xlab(side_long) +
            ylab('Base type') +
            scale_fill_distiller(palette = 'PuOr', name = 'log10(probability of deletion)', limits = limits, na.value = 'gray80') + 
            theme(text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 18), axis.text.x = element_blank(), legend.position = 'bottom', legend.direction = 'horizontal', legend.justification="center", axis.title = element_text(size = 24)) +
            guides(fill = guide_colourbar(barwidth = 35, barheight = 2)) 
        
        assign(name, plot)
    }
    return(list('5end_plot' = get('5endplot'), '3end_plot' = get('3endplot')))
}




v_motif_heatmap = plot_oriented_motif_coefficient_heatmap(coefs[trim_type == 'v_trim'], orient='5to3', limits = bounds)
j_motif_heatmap = plot_oriented_motif_coefficient_heatmap(coefs[trim_type == 'j_trim'], orient='3to5', limits = bounds)

v_motif_heatmap = v_motif_heatmap + theme(legend.position = 'none') 
j_motif_heatmap = j_motif_heatmap + theme(legend.position = 'none') 

v_base_plots = plot_indiv_base_count_coefficient(coefs[trim_type == 'v_trim'], limits=bounds)
j_base_plots = plot_indiv_base_count_coefficient(coefs[trim_type == 'j_trim'], limits=bounds)

v_5end_base_plot = v_base_plots[['5end_plot']] + theme(legend.position = 'none') 
v_3end_base_plot = v_base_plots[['3end_plot']] + theme(legend.position = 'none') 
j_5end_base_plot = j_base_plots[['5end_plot']] + theme(legend.position = 'none') 
j_3end_base_plot = j_base_plots[['3end_plot']] + theme(legend.position = 'none') 

mh_heatmap = plot_ligation_mh_coefficient_heatmap_single_group(coefs, with_values = FALSE, limits = bounds) + xlab('Number of MH nucleotides\nin ligation configuration') + theme(text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_text(size = 24), axis.title.y = element_blank())

config_heatmap = plot_average_mh_coefficient_heatmap_single_group(coefs, with_values = FALSE, limits = bounds) + xlab('Average number of MH nucleotides\nacross possible ligation configurations') + theme(text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_text(size = 24), axis.title.y = element_blank())

config_heatmap = config_heatmap + theme(legend.position = 'none')

# isolate legend
mh_heatmap = mh_heatmap + guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = 38, barheight = 2)) 
legend = get_legend(mh_heatmap)
mh_heatmap = mh_heatmap + theme(legend.position = 'none')

ggsave(paste0(file_root, '/vmotif.pdf'), plot = v_motif_heatmap, width = 6, height = 5.5, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)
ggsave(paste0(file_root, '/jmotif.pdf'), plot = j_motif_heatmap, width = 6, height = 5.5, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)
ggsave(paste0(file_root, '/v5_base_count.pdf'), plot = v_5end_base_plot, width = 4.5, height = 5.5, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)
ggsave(paste0(file_root, '/v3_base_count.pdf'), plot = v_3end_base_plot, width = 4.5, height = 5.5, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)
ggsave(paste0(file_root, '/j5_base_count.pdf'), plot = j_5end_base_plot, width = 4.5, height = 5.5, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)
ggsave(paste0(file_root, '/j3_base_count.pdf'), plot = j_3end_base_plot, width = 4.5, height = 5.5, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)
ggsave(paste0(file_root, '/trim_mh.pdf'), plot = config_heatmap, width = 7, height = 3, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)
ggsave(paste0(file_root, '/lig_mh.pdf'), plot = mh_heatmap, width = 7, height = 3, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)
ggsave(paste0(file_root, '/legend.pdf'), plot = legend, width = 8, height = 2, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE
)
