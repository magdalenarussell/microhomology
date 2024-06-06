source('config/config.R')
source('config/file_paths.R')

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(Biostrings)
library(mclogit)
library(matrixcalc)
library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)

args = commandArgs(trailingOnly=TRUE)

ANNOTATION_TYPE <<- 'igor_mh_sim_alpha'
PARAM_GROUP <<- args[1]
stopifnot(PARAM_GROUP %in% c('nonproductive_v-j_trim_ligation-mh', 'both_v-j_trim_ligation-mh'))
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))
SAMPLE_ANNOT <<- TRUE
NCPU <<- as.numeric(args[2])
# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 1
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 2
MODEL_TYPE <<- 'motif_two-side-base-count-beyond_mh-config-count_ligation-mh'

TRIMMING_PROB_TYPE <<- 'uniform'
stopifnot(TRIMMING_PROB_TYPE %in% c('igor', 'uniform'))

LIGATION_MH_PARAM <<- 0
INT_MH_PARAM <<- 0

NP_COND <<- args[3]

if (PARAM_GROUP %like% 'both'){
    stopifnot(NP_COND == FALSE)
}

source(paste0(MOD_PROJECT_PATH,'/analysis_scripts/trim_lig_config_mh_simulator/ligation-mh_simulator_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

# total number of sequences
total = 2500000

# get V and J gene choice probs
gene_choice_probs = get_igor_gene_usage_params()
whole_nucseq = get_oriented_whole_nucseqs()
gene_choice_probs$v_choice = merge(gene_choice_probs$v_choice, whole_nucseq[, c('gene', 'v_gene_sequence')], by.x = 'v_gene', by.y = 'gene') 
gene_choice_probs$j_choice = merge(gene_choice_probs$j_choice, whole_nucseq[, c('gene', 'j_gene_sequence')], by.x = 'j_gene', by.y = 'gene') 

# get all possible configurations
configs = read_frames_data()
configs = configs[v_trim <= UPPER_TRIM_BOUND & j_trim <= UPPER_TRIM_BOUND]

if (NP_COND == TRUE){
    configs = configs[frame_type == 'Out' | frame_stop == TRUE]
}

set.seed(123) 

registerDoParallel(cores=NCPU)
sim_data = foreach(subset = seq(total/100), .combine=rbind) %dopar% {
    subset_sim = data.table()
    for (i in seq(100)){
        # print(paste0('starting sim ', i))
        # sample a V and J gene
        temp = cbind(gene_choice_probs$v_choice[sample(.N, 1, prob = v_gene_prob)],
                     gene_choice_probs$j_choice[sample(.N, 1, prob = j_gene_prob)]) 

        temp_configs = configs[v_gene == temp$v_gene & j_gene == temp$j_gene]

        # sample a V and J gene trimming amount
        trim_cols = c('v_gene', 'j_gene', 'v_trim', 'j_trim')
        temp_trims = unique(temp_configs[, ..trim_cols])

        if (nrow(temp_trims) == 0){
            next
        }

        cols = c('v_gene', 'v_gene_prob', 'v_gene_sequence', 'j_gene', 'j_gene_prob', 'j_gene_sequence')
        temp = temp[, ..cols]

        ## Create a probability vector where lower values have higher probability
        ## For example, using linearly decreasing weights
        temp_trim = temp_trims[sample(.N, 1)]

        temp = merge(temp, temp_trim, by = c('v_gene', 'j_gene'))

        temp_configs2 = temp_configs[v_gene == temp$v_gene & j_gene == temp$j_gene & v_trim == temp$v_trim & j_trim == temp$j_trim]

        lig_options = unique(temp_configs2$ligation_mh)
        if (length(lig_options) > 1){
            lig_draw = sample(lig_options, 1)
        } else {
            lig_draw = lig_options
        }

        temp$ligation_mh = as.numeric(lig_draw)
        subset_sim = rbind(subset_sim, temp)
    }
    subset_sim
}
stopImplicitCluster()

cols = colnames(sim_data)[colnames(sim_data) != 'ligation_attempt']

condensed_sim = sim_data[, .N, by = cols]
setnames(condensed_sim, 'N', 'count')

condensed_sim$vj_insert = 0

motif_data = get_all_nuc_contexts(condensed_sim, 'mh_sim', gene_type = GENE_NAME, trim_type = TRIM_TYPE)
ncount = sum(motif_data$count)

# fill in missing sites
motif_data = merge(motif_data, configs)
stopifnot(sum(motif_data$count) == ncount)

filled_motif_data = filter_motif_data_for_possible_sites(motif_data)

# filter for productivity
if (grepl('nonprod', PARAM_GROUP)){
    filled_motif_data[frame_type == 'In' & frame_stop == FALSE, count := 0]
}

processed0 = inner_aggregation_processing(filled_motif_data, gene_type = GENE_NAME, trim_type = TRIM_TYPE, only_nonprod_sites = FALSE)
processed = subset_processed_data(processed0, trim_type = TRIM_TYPE, gene_type = GENE_NAME)

old_ANNOTATION_TYPE <<- ANNOTATION_TYPE
ANNOTATION_TYPE <<- paste0(old_ANNOTATION_TYPE, '_from_', TRIMMING_PROB_TYPE, '_MHprob', LIGATION_MH_PARAM, '_trimMHprob', INT_MH_PARAM, '_simple_NPconditioned', NP_COND)

filename = processed_data_path(sample_annotation = SAMPLE_ANNOT)
fwrite(processed, filename, sep = '\t')

if (grepl('nonprod', PARAM_GROUP)){
    processed1 = inner_aggregation_processing(filled_motif_data, gene_type = GENE_NAME, trim_type = TRIM_TYPE, only_nonprod_sites = TRUE)
    processed = subset_processed_data(processed1, trim_type = TRIM_TYPE, gene_type = GENE_NAME)

    filename = processed_data_path(sample_annotation = SAMPLE_ANNOT)
    fwrite(processed, filename, sep = '\t')
}

if (length(colnames(processed)[colnames(processed) %like% 'ligation']) > 1){
    processed = processed[, -25]
}

processed$vj_insert = 0
# now get all annotations
all = get_all_annotations(processed)

# Get oriented full sequences and group genes by common features
all = get_oriented_full_sequences(all, gene_type = GENE_NAME)

all$vj_insert = 0
cols = c('v_gene', 'j_gene', 'v_trim', 'j_trim', 'ligation_mh', 'vj_insert', 'v_gene_sequence', 'j_gene_sequence', 'processed_sequence')
condensed_sim = all[, sum(count), by = cols]
setnames(condensed_sim, 'V1', 'count')
condensed_sim[, index := .GRP, by = .(v_gene, j_gene, processed_sequence)]

motif_data = get_all_nuc_contexts(condensed_sim, 'adjusted-mh_sim', gene_type = GENE_NAME, trim_type = TRIM_TYPE)

# reset annotation type
ANNOTATION_TYPE <<- 'igor_mh_sim_alpha'
SAMPLE_ANNOT <<- FALSE
source(paste0(MOD_PROJECT_PATH,'/analysis_scripts/ligation-mh_signal_simulator/ligation-mh_simulator_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

# fill in missing sites
filled_motif_data = filter_motif_data_for_possible_sites(motif_data)

processed = inner_aggregation_processing(filled_motif_data, gene_type = GENE_NAME, trim_type = TRIM_TYPE, only_nonprod_sites = FALSE)
processed = subset_processed_data(processed, trim_type = TRIM_TYPE, gene_type = GENE_NAME)

old_ANNOTATION_TYPE <<- ANNOTATION_TYPE
ANNOTATION_TYPE <<- paste0(old_ANNOTATION_TYPE, '_from_', TRIMMING_PROB_TYPE, '_MHprob', LIGATION_MH_PARAM, '_trimMHprob', INT_MH_PARAM, '_simple_NPconditioned', NP_COND)

filename = processed_data_path(sample_annotation = SAMPLE_ANNOT)
fwrite(processed, filename, sep = '\t')

if (grepl('nonprod', PARAM_GROUP)){
    processed1 = inner_aggregation_processing(filled_motif_data, gene_type = GENE_NAME, trim_type = TRIM_TYPE, only_nonprod_sites = TRUE)
    processed = subset_processed_data(processed1, trim_type = TRIM_TYPE, gene_type = GENE_NAME)

    filename = processed_data_path(sample_annotation = SAMPLE_ANNOT)
    fwrite(processed, filename, sep = '\t')
}
