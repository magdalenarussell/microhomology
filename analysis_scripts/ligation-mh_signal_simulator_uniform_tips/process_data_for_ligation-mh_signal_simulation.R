source('mechanistic-trimming/config/config.R')
source('mechanistic-trimming/config/file_paths.R')

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
PARAM_GROUP <<- 'nonproductive_v-j_trim_ligation-mh'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))
SAMPLE_ANNOT <<- TRUE
NCPU <<- as.numeric(args[1])
# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 1
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 2
MODEL_TYPE <<- 'motif_two-side-base-count-beyond_ligation-mh'

TRIMMING_PROB_TYPE <<- 'uniform_tips'

LIGATION_MH_PARAM <<- as.numeric(args[2])

source(paste0(MOD_PROJECT_PATH,'/analysis_scripts/ligation-mh_signal_simulator_scripts/ligation-mh_simulator_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

# total number of sequences
total = 2500000

# get V and J gene choice probs
gene_choice_probs = get_igor_gene_usage_params()
whole_nucseq = get_oriented_whole_nucseqs()
gene_choice_probs$v_choice = merge(gene_choice_probs$v_choice, whole_nucseq[, c('gene', 'v_gene_sequence')], by.x = 'v_gene', by.y = 'gene') 
gene_choice_probs$j_choice = merge(gene_choice_probs$j_choice, whole_nucseq[, c('gene', 'j_gene_sequence')], by.x = 'j_gene', by.y = 'gene') 

# get all ligation probabilities
all_lig_probs = get_ligation_probabilities(LIGATION_MH_PARAM, seq(0, 15))
all_lig_probs = all_lig_probs[names(all_lig_probs) != 'no_ligation']

# get all possible configurations
configs = read_frames_data()

set.seed(123) 

registerDoParallel(cores=NCPU)
sim_data = foreach(subset = seq(total/100), .combine=rbind) %dopar% {
    subset_sim = data.table()
    for (i in seq(100)){
        # print(paste0('starting sim ', i))
        # sample a V and J gene
        temp = cbind(gene_choice_probs$v_choice[sample(.N, 1, prob = v_gene_prob)],
                     gene_choice_probs$j_choice[sample(.N, 1, prob = j_gene_prob)]) 

        # sample a V and J gene trimming amount
        cols = c('v_gene', 'j_gene', 'v_trim', 'j_trim', 'ligation_mh')
        possible_trims = configs[v_gene == temp$v_gene & j_gene == temp$j_gene][, ..cols]
        possible_trims$trim_prob = 1

        if (nrow(possible_trims) == 0){
            next
        }

        cols = c('v_gene', 'v_gene_prob', 'v_gene_sequence', 'j_gene', 'j_gene_prob', 'j_gene_sequence')
        temp = temp[, ..cols]

        ## Create a probability vector where lower values have higher probability
        ## For example, using linearly decreasing weights
        temp_trim = possible_trims[sample(.N, 1, prob=trim_prob)]

        temp = merge(temp, temp_trim, by = c('v_gene', 'j_gene'))

        temp_configs = configs[v_gene == temp$v_gene & j_gene == temp$j_gene & v_trim == temp$v_trim & j_trim == temp$j_trim]

        lig_probs = all_lig_probs[names(all_lig_probs) %in% unique(temp_configs$ligation_mh)]

        lig_draw = sample(names(lig_probs), 1, prob = lig_probs)

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

# Because there is no selection effect in these simulated data, I am not
# restricting to nonproductive sequences! 
# fill in missing sites
cols2 = c('v_gene', 'j_gene', 'v_trim', 'j_trim', 'ligation_mh')
tog = merge(motif_data, configs, by = cols2)

filled_motif_data = fill_in_missing_possible_sites(configs, tog, trim_type = TRIM_TYPE, gene_type = GENE_NAME)

processed = inner_aggregation_processing(filled_motif_data, gene_type = GENE_NAME, trim_type = TRIM_TYPE)
processed = subset_processed_data(processed, trim_type = TRIM_TYPE, gene_type = GENE_NAME)

old_ANNOTATION_TYPE <<- ANNOTATION_TYPE
ANNOTATION_TYPE <<- paste0(old_ANNOTATION_TYPE, '_from_', TRIMMING_PROB_TYPE, '_MHprob', LIGATION_MH_PARAM)

filename = processed_data_path(sample_annotation = SAMPLE_ANNOT)
fwrite(processed, filename, sep = '\t')
