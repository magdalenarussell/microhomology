source('config/config.R')

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(Biostrings)
library(mclogit)
library(matrixcalc)
library(RhpcBLASctl)
library(cowplot)
omp_set_num_threads(1)
blas_set_num_threads(1)

args = commandArgs(trailingOnly=TRUE)

ANNOTATION_TYPE <<- 'igor'
TRIM_TYPE <<- 'v_trim'
PRODUCTIVITY <<- 'nonproductive'
MOTIF_TYPE <<- 'unbounded' 
NCPU <<- 2
GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
MODEL_GROUP <<- 'all_subjects'
GENE_WEIGHT_TYPE <<- 'p_gene_marginal' 
# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 1
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 2
UPPER_TRIM_BOUND <<- 14
LOWER_TRIM_BOUND <<- 2
MODEL_TYPE <<- 'motif_two-side-base-count-beyond'
LEFT_SIDE_TERMINAL_MELT_LENGTH <<- 10

# load BCR data
VALIDATION_DATA_DIR <<- args[1]
VALIDATION_TYPE <<- 'validation_data_igh'
VALIDATION_TRIM_TYPE <<- args[2]
VALIDATION_PRODUCTIVITY <<- 'nonproductive'
VALIDATION_GENE_NAME <<- paste0(substring(VALIDATION_TRIM_TYPE, 1, 1), '_gene')

source(paste0(MOD_PROJECT_PATH,'scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'scripts/model_fitting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'plotting_scripts/plotting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'plotting_scripts/residual_comparison_functions.R'))
source(paste0(MOD_PROJECT_PATH,'scripts/analysis_scripts/pwm_profile_functions.R'))

# load model
if (MODEL_TYPE != 'null') {
    model = load_model()
} else {
    model = 'null'
} 

# get model coefficient data
pwm = get_model_coefficient_data()

# get predicted trimming distributions for training data, and calculate subject-gene weight
source(paste0('scripts/sampling_procedure_functions/p_gene_given_subject.R'), local = TRUE)
predicted_trims = get_predicted_distribution_data() 
predicted_trims = calculate_subject_gene_weight(predicted_trims)

# get PWM score for each observation
training_pwm_scores = get_all_pwm_score(predicted_trims, pwm, NULL)
training_avg = compare_weight_by_motif(predicted_trims)
training = merge(training_avg, unique(training_pwm_scores[, -c('trim_length')]), by = c('motif', 'gene'))

# set parameters for IGH data set
ANNOTATION_TYPE <<- VALIDATION_TYPE
TYPE <<- 'validation_data' 
TRIM_TYPE <<- VALIDATION_TRIM_TYPE
GENE_NAME <<- VALIDATION_GENE_NAME
PRODUCTIVITY <<- VALIDATION_PRODUCTIVITY

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('scripts/model_evaluation_functions.R')

# get IGH data and predict distributions using model
validation_data = aggregate_validation_data(directory = VALIDATION_DATA_DIR)
validation_data$predicted_prob = temp_predict(model, newdata = validation_data)
validation_data[, empirical_prob := count/sum(count), by = .(subject, gene)]

# calculate subject-gene weight
source(paste0('scripts/sampling_procedure_functions/p_gene_given_subject.R'), local = TRUE)
validation_data = calculate_subject_gene_weight(validation_data)

# get PWM score for each observation
validation_pwm_scores = get_all_pwm_score(validation_data, pwm, NULL)
validation_avg = compare_weight_by_motif(validation_data)
validation = merge(validation_avg, unique(validation_pwm_scores[, -c('trim_length')]), by = c('motif', 'gene'))

#get plot path
ANNOTATION_TYPE <<- 'igor'
path = get_manuscript_path()
path = file.path(path, 'motif_analysis')
dir.create(path)

# combine data sets
together = rbind(training, validation)
together[gene %like% 'IGH', type := 'IGH']
together[gene %like% 'TRB', type := 'TRB']

# get motif frequencies
freqs = together[, .N, by = .(motif, type, pwm_score)]
freqs[, total := sum(N), by = type]
freqs[, freq := N/total]

# order frequencies
ordered_freqs = unique(freqs[order(freqs$pwm_score)]$motif) 

# create an empty data set for spacing purposes
spacing = data.table(motif = c('', ''), type = c('TRB', 'IGH'), pwm_score = c(NA, NA), freq = c(0.09, 0.09))

# combine data sets, re-orient IGH data for plotting purposes
freqs2 = rbind(freqs, spacing, fill = TRUE)
freqs2$motif = factor(freqs2$motif, levels = ordered_freqs)
freqs2[is.na(pwm_score), motif := '']
freqs2[type == 'IGH', freq := -1*freq]

# plot frequencies
plot = ggplot(freqs2, aes(y = freq, x = motif, fill = pwm_score)) +
    geom_bar(position = position_dodge(width=1), stat='identity') + 
    facet_share(~type, dir = "h", scales = "free", reverse_num = TRUE) +
    coord_flip() +
    xlab('') +
    ylab('Germline motif frequency') +
    scale_fill_distiller(palette = 'PuOr', name = 'PWM motif score', limits = c(-1.24, 1.24), na.value="transparent") +
    theme_cowplot(font_family = 'Arial') + 
    theme(text = element_text(size = 25), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5), axis.text.x = element_text(size = 20)) + 
    guides(fill = guide_colourbar(barwidth = 2, barheight = 15))+
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) 

# save plot
ggsave(file.path(path, 'germline_motif_freq.pdf'), plot = plot, width = 16, height = 45, units = 'in', dpi = 750, device = cairo_pdf)


