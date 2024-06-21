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
PARAM_GROUP <<- 'nonproductive_v-j_trim'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))
NCPU <<- as.numeric(args[1])
# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 1
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 2
MODEL_TYPE <<- 'motif_two-side-base-count-beyond_interior-mh-count'

INT_MH_PARAM <<- as.numeric(args[2])

source(paste0(MOD_PROJECT_PATH,'/mh_simulation_scripts/intermediate-mh_signal_simulator_scripts/intermediate-mh_simulator_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

stopifnot(LOCUS == 'alpha')

# total number of sequences
total = 4000000

# get V and J gene choice probs
gene_choice_probs = get_igor_gene_usage_params()
whole_nucseq = get_oriented_whole_nucseqs()
gene_choice_probs$v_choice = merge(gene_choice_probs$v_choice, whole_nucseq[, c('gene', 'v_gene_sequence')], by.x = 'v_gene', by.y = 'gene') 
gene_choice_probs$j_choice = merge(gene_choice_probs$j_choice, whole_nucseq[, c('gene', 'j_gene_sequence')], by.x = 'j_gene', by.y = 'gene') 

# get Vtrim and Jtrim probs
trim_probs = get_trimming_probs()

# get all possible configurations
configs = get_all_possible_configs(gene_choice_probs, trim_probs, int_mh_prob=INT_MH_PARAM)

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
        possible_trims = configs[v_gene == temp$v_gene & j_gene == temp$j_gene]

        if (nrow(possible_trims) == 0){
            next
        }

        sampled = possible_trims[sample(.N, 1, prob=mh_adj_joint_trim_prob)]
        subset_sim = rbind(subset_sim, sampled)
    }
    subset_sim
}
stopImplicitCluster()

cols = colnames(sim_data)
condensed_sim = sim_data[, .N, by = cols]
setnames(condensed_sim, 'N', 'count')

condensed_sim$vj_insert = 5

motif_data = get_all_nuc_contexts(condensed_sim, 'mh_sim', gene_type = GENE_NAME, trim_type = TRIM_TYPE)

processed = inner_aggregation_processing(motif_data, gene_type = GENE_NAME, trim_type = TRIM_TYPE)
processed = subset_processed_data(processed, trim_type = TRIM_TYPE, gene_type = GENE_NAME)

old_ANNOTATION_TYPE <<- ANNOTATION_TYPE
ANNOTATION_TYPE <<- paste0(old_ANNOTATION_TYPE, '_from_uniform_MHprob', INT_MH_PARAM)

filename = processed_data_path(sample_annotation = SAMPLE_ANNOT)
fwrite(processed, filename, sep = '\t')
