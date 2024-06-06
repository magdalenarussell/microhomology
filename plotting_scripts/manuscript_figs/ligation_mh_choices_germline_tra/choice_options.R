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

MODEL_TYPE <<- 'motif_two-side-base-count-beyond_mh-config-count_ligation-mh'

L2 <<- 'False'

source(paste0(MOD_PROJECT_PATH,'/mh_simulation_scripts/ligation-mh_signal_simulator/ligation-mh_simulator_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

# get all possible configurations
configs = read_frames_data()
configs = configs[v_trim <= UPPER_TRIM_BOUND & j_trim <= UPPER_TRIM_BOUND]
gene_choice_probs = get_igor_gene_usage_params()
top_V = gene_choice_probs$v_choice[order(-v_gene_prob)][1:6]
top_J = gene_choice_probs$j_choice[order(-j_gene_prob)][1:6]

all_options_df = configs[v_gene %in% top_V$v_gene & j_gene %in% top_J$j_gene, .N, by = .(v_gene, j_gene, v_trim, j_trim)]

zero_mh_options = configs[ligation_mh == 0][v_gene %in% top_V$v_gene & j_gene %in% top_J$j_gene]

all_options_df = merge(all_options_df, zero_mh_options[, c('v_gene', 'j_gene', 'v_trim', 'j_trim', 'ligation_mh')], by = c('v_gene', 'j_gene', 'v_trim', 'j_trim'), all.x = TRUE)

all_options_plot = ggplot(all_options_df) +
                   facet_grid(cols = vars(j_gene), rows = vars(v_gene))+
                   geom_tile(aes(x = v_trim, y = j_trim, fill = N)) +
                   geom_point(data = all_options_df[ligation_mh == 0], aes(x = v_trim, y = j_trim, color = ligation_mh)) +
                   ylab("J-gene trimming length") +
                   xlab("V-gene trimming length") +
                   theme_cowplot(font_family = 'Arial') + 
                   background_grid(major = 'xy') + 
                   panel_border(color = 'gray60', size = 1.5) +
                   scale_fill_viridis_c(option = 'rocket', direction = -1, limits = c(1, 4))+
                   theme(text = element_text(size = 40), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 30), legend.key.size = unit(3, "cm"))+
                   labs(fill = 'Number of ligation\nMH choices', color = 'Zero MH\nligation option')



# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/ligation_mh_choices_germline_tra/all_choices_heatmap_group_top_genes.pdf')

ggsave(file_name, plot = all_options_plot, width = 35, height = 30, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE) 

# save data               
fwrite(all_options_df, paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/ligation_mh_choices_germline_tra/all_choices_top_genes.tsv'), sep = '\t') 

# get possible configs after filtering for nonprods
configs_nonprod = configs[frame_type == 'Out' | frame_stop == TRUE]
zero_mh_options = configs_nonprod[ligation_mh == 0][v_gene %in% top_V$v_gene & j_gene %in% top_J$j_gene]
all_options_df = configs_nonprod[v_gene %in% top_V$v_gene & j_gene %in% top_J$j_gene, .N, by = .(v_gene, j_gene, v_trim, j_trim)]

all_options_df = merge(all_options_df, zero_mh_options[, c('v_gene', 'j_gene', 'v_trim', 'j_trim', 'ligation_mh')], by = c('v_gene', 'j_gene', 'v_trim', 'j_trim'), all.x = TRUE)

all_options_plot = ggplot(all_options_df) +
                   facet_grid(cols = vars(j_gene), rows = vars(v_gene))+
                   geom_tile(aes(x = v_trim, y = j_trim, fill = N)) +
                   geom_point(data = all_options_df[ligation_mh == 0], aes(x = v_trim, y = j_trim, color = ligation_mh)) +
                   ylab("J-gene trimming length") +
                   xlab("V-gene trimming length") +
                   theme_cowplot(font_family = 'Arial') + 
                   background_grid(major = 'xy') + 
                   panel_border(color = 'gray60', size = 1.5) +
                   scale_fill_viridis_c(option = 'rocket', direction = -1, limits=c(1, 4))+
                   theme(text = element_text(size = 40), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 30), legend.key.size = unit(3, "cm"))+
                   labs(fill = 'Number of ligation\nMH choices', color = 'Zero MH\nligation option')


# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/ligation_mh_choices_germline_tra/all_choices_heatmap_group_top_genes_nonprod.pdf')

ggsave(file_name, plot = all_options_plot, width = 35, height = 30, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE) 

# save data               
fwrite(all_options_df, paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/ligation_mh_choices_germline_tra/all_choices_top_genes_nonprod.tsv'), sep = '\t') 
