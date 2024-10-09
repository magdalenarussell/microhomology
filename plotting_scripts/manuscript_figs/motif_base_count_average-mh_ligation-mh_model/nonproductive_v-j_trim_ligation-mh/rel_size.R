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
pval_path = get_bootstrap_pvalue_path(L2)
coefs = fread(pval_path)
coefs$bonf_level = 0.05/nrow(coefs)
coefs[pvalue < bonf_level, significant := TRUE]
coefs[pvalue >= bonf_level, significant := FALSE]
coefs[value > 0, star_position := value/log(10)-0.02]
coefs[value < 0, star_position := value/log(10)- 0.05]

# make clean names
coefs[coefficient == 'ligation_mh', clean_name := 'ligation-related MH']
coefs[coefficient == 'average_mh', clean_name := 'trimming-related MH']
coefs[coefficient == 'base_count' & side %like% '5'& trim_type %like% 'v', clean_name := paste('5\' V-gene ', base, ' count')]
coefs[coefficient == 'base_count' & side %like% '5'& trim_type %like% 'j', clean_name := paste('5\' J-gene ', base, ' count')]
coefs[coefficient == 'base_count' & side %like% '3'& trim_type %like% 'v', clean_name := paste('3\' V-gene ', base, ' count')]
coefs[coefficient == 'base_count' & side %like% '3'& trim_type %like% 'j', clean_name := paste('3\' J-gene ', base, ' count')]
coefs[coefficient == 'motif' & side %like% '5'& trim_type %like% 'v', clean_name := paste('5\' V-motif ', base)]
coefs[coefficient == 'motif' & side %like% '5'& trim_type %like% 'j', clean_name := paste('5\' J-motif ', base)]
coefs[coefficient == 'motif' & side %like% '3'& trim_type %like% 'v', clean_name := paste('3\'-', position, ' V-motif ', base)]
coefs[coefficient == 'motif' & side %like% '3'& trim_type %like% 'j', clean_name := paste('3\'-', position, ' J-motif ', base)]

name_order = coefs[order(-abs(value))]$clean_name
coefs$clean_name = factor(coefs$clean_name, levels = name_order)

coefs[coefficient %like% 'mh', type := 'MH']
coefs[coefficient %like% 'base' & side %like% '5', type := '5\' base count']
coefs[coefficient %like% 'base' & side %like% '3', type := '3\' base count']
coefs[coefficient %like% 'motif', type := 'Motif']

colormap = c('5\' base count'='#66a61e', 'Motif' = '#e7298a', '3\' base count'='#666666', 'MH'='#e6ab02')

coefs$type = factor(coefs$type, levels = names(colormap))

bar = ggplot(coefs, aes(x = clean_name, y = value/log(10), fill = type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colormap, name = 'Parameter type') +
  # scale_fill_distiller(palette = 'PuOr', name = 'log10(probability of deletion)', limits = c(-0.450, 0.450), na.value = 'gray80') + 
  scale_y_continuous(expand = c(0, 0), limits = c(min(coefs$value/log(10)) - 0.05, max(coefs$value/log(10)) + 0.05)) +
  geom_text(aes(label = ifelse(significant, "*", ""), y = star_position), vjust = 0, size = 3, color = "black") +
  theme_cowplot(font_family = 'Arial') + 
  xlab('Parameter') +
  ylab(expression(log[10]("parameter value"))) +
  background_grid(major = 'xy') + 
  panel_border(color = 'gray60', size = 1.5) +
  theme(text = element_text(size = 12), axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(size = 10), axis.text.x = element_text(angle = -45, hjust = 0, size = 10), axis.title = element_text(size = 12), plot.margin = unit(c(0.5,2.2,0.5,0.5), "cm"), legend.position = 'none') 

ggsave(paste0(file_root, '/rel_effect.pdf'), plot = bar, width = 10.2, height = 3, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE)
