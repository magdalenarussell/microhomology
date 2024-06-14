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
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))
SAMPLE_ANNOT <<- TRUE
NCPU <<- as.numeric(args[2])
# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 1
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 2
MODEL_TYPE <<- 'motif_two-side-base-count-beyond_mh-config-count_ligation-mh'

TRIMMING_PROB_TYPE <<- args[3]
stopifnot(TRIMMING_PROB_TYPE %in% c('igor', 'motif_two-side-base-count-beyond', 'uniform', 'mh_adjusted_motif-two-side-base-count-beyond'))

LIGATION_MH_PARAM <<- as.numeric(args[4])
INT_MH_PARAM <<- as.numeric(args[5])

NP_COND <<- args[6]

if (PARAM_GROUP %like% 'both'){
    stopifnot(NP_COND == FALSE)
    params = c('both_v-j_trim_ligation-mh', 'nonproductive_v-j_trim_ligation-mh')
} else {
    params = PARAM_GROUP
}

source(paste0(MOD_PROJECT_PATH,'/mh_simulation_scripts/ligation-mh_signal_simulator/ligation-mh_simulator_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

stopifnot(LOCUS == 'alpha')

# total number of sequences
total = 2500000

# get V and J gene choice probs
gene_choice_probs = get_igor_gene_usage_params()
whole_nucseq = get_oriented_whole_nucseqs()
gene_choice_probs$v_choice = merge(gene_choice_probs$v_choice, whole_nucseq[, c('gene', 'v_gene_sequence')], by.x = 'v_gene', by.y = 'gene') 
gene_choice_probs$j_choice = merge(gene_choice_probs$j_choice, whole_nucseq[, c('gene', 'j_gene_sequence')], by.x = 'j_gene', by.y = 'gene') 

# get all possible configurations
configs = read_frames_data()
if (NP_COND == TRUE){
    configs = configs[frame_type == 'Out' | frame_stop == TRUE]
    file_name = frame_data_path()
    if (!file.exists(file_name)){
        file_name = str_replace(file_name, 'TRA', 'TRA_NPcond')
        dir.create(dirname(file_name), recursive = TRUE)
        fwrite(configs, file_name, sep = '\t')
    }
}

# get Vtrim and Jtrim probs
trim_probs = get_trimming_probs(TRIMMING_PROB_TYPE, INT_MH_PARAM, configs)

# get all ligation probabilities
all_lig_probs = get_ligation_probabilities(LIGATION_MH_PARAM, seq(0, 15))

set.seed(1) 
registerDoParallel(cores=NCPU)
sim_data = foreach(subset = seq(total/100), .combine=rbind) %dopar% {
    subset_sim = data.table()
    for (i in seq(100)){
        # sample a V and J gene
        temp = cbind(gene_choice_probs$v_choice[sample(.N, 1, prob = v_gene_prob)],
                     gene_choice_probs$j_choice[sample(.N, 1, prob = j_gene_prob)]) 
        cols = c('v_gene', 'v_gene_prob', 'v_gene_sequence', 'j_gene', 'j_gene_prob', 'j_gene_sequence')
        temp = temp[, ..cols]

        # get possible trimming configurations
        possible_trims = trim_probs[v_gene == temp$v_gene & j_gene == temp$j_gene]
        if (nrow(possible_trims) == 0){
            next
        }
        if (all(unique(possible_trims$mh_adj_joint_trim_prob) == 0)){
            next
        }

        # sample a trimming configuration according to MH-adjusted trimming probability
        temp_trim = possible_trims[sample(.N, 1, prob=mh_adj_joint_trim_prob)]
        temp = merge(temp, temp_trim, by = c('v_gene', 'j_gene'))

        # get all possible ligation configurations
        temp_configs = configs[v_gene == temp$v_gene & j_gene == temp$j_gene & v_trim == temp$v_trim & j_trim == temp$j_trim]

        # sample ligation configuration according to MH-adjusted ligation probability
        lig_probs = all_lig_probs[names(all_lig_probs) %in% unique(temp_configs$ligation_mh)]
        lig_options = names(lig_probs)
        if (length(lig_options) > 1){
            lig_draw = sample(lig_options, 1, prob = lig_probs)
        } else {
            lig_draw = lig_options
        }

        temp$ligation_mh = as.numeric(lig_draw)
        subset_sim = rbind(subset_sim, temp)
    }
    subset_sim
}
stopImplicitCluster()

for (PARAM_GROUP in params){
    ANNOTATION_TYPE <<- 'igor_mh_sim_alpha'
    
    source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))
    SAMPLE_ANNOT <<- TRUE
    source(paste0(MOD_PROJECT_PATH,'/mh_simulation_scripts/ligation-mh_signal_simulator/ligation-mh_simulator_functions.R'))
    source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

    ###############################################
    ### process data using observed simulations ###
    ###############################################

    cols = colnames(sim_data)[colnames(sim_data) != 'ligation_attempt']

    # condense data
    condensed_sim = sim_data[, .N, by = cols]
    setnames(condensed_sim, 'N', 'count')
    condensed_sim$vj_insert = 0

    # get nucleotide context
    motif_data = get_all_nuc_contexts(condensed_sim, 'mh_sim', gene_type = GENE_NAME, trim_type = TRIM_TYPE)
    ncount = sum(motif_data$count)

    # fill in missing sites
    motif_data = merge(motif_data, configs)
    stopifnot(sum(motif_data$count) == ncount)
    filled_motif_data = filter_motif_data_for_possible_sites(motif_data)

    # get model parameters and subset data
    processed0 = inner_aggregation_processing(filled_motif_data, gene_type = GENE_NAME, trim_type = TRIM_TYPE, only_nonprod_sites = ONLY_NONPROD_SITES)
    processed = subset_processed_data(processed0, trim_type = TRIM_TYPE, gene_type = GENE_NAME)

    # save data
    old_ANNOTATION_TYPE <<- ANNOTATION_TYPE
    ANNOTATION_TYPE <<- paste0(old_ANNOTATION_TYPE, '_from_', TRIMMING_PROB_TYPE, '_MHprob', LIGATION_MH_PARAM, '_trimMHprob', INT_MH_PARAM)
    if (NP_COND == TRUE){
        ANNOTATION_TYPE <<- paste0(ANNOTATION_TYPE, '_NPcond')
    }

    filename = processed_data_path(sample_annotation = SAMPLE_ANNOT)
    fwrite(processed, filename, sep = '\t')

    ######################################################
    ### get all annotations using observed simulations ###
    ######################################################

    # now get all annotations
    processed$vj_insert = 0
    all = get_all_annotations(processed)

    # Get oriented full sequences and group genes by common features
    all = get_oriented_full_sequences(all, gene_type = GENE_NAME)

    # condense data
    all$vj_insert = 0
    cols = c('v_gene', 'j_gene', 'v_trim', 'j_trim', 'ligation_mh', 'vj_insert', 'v_gene_sequence', 'j_gene_sequence', 'processed_sequence')
    condensed_sim = all[, sum(count), by = cols]
    setnames(condensed_sim, 'V1', 'count')
    condensed_sim[, index := .GRP, by = .(v_gene, j_gene, processed_sequence)]

    # get nuc context
    motif_data = get_all_nuc_contexts(condensed_sim, 'adjusted-mh_sim', gene_type = GENE_NAME, trim_type = TRIM_TYPE)

    # reset annotation type
    ANNOTATION_TYPE <<- 'igor_mh_sim_alpha'
    SAMPLE_ANNOT <<- FALSE
    source(paste0(MOD_PROJECT_PATH,'/mh_simulation_scripts/ligation-mh_signal_simulator/ligation-mh_simulator_functions.R'))
    source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

    # fill in missing sites
    filled_motif_data = filter_motif_data_for_possible_sites(motif_data)

    # get model parameters and subset data
    processed = inner_aggregation_processing(filled_motif_data, gene_type = GENE_NAME, trim_type = TRIM_TYPE, only_nonprod_sites = ONLY_NONPROD_SITES)
    processed = subset_processed_data(processed, trim_type = TRIM_TYPE, gene_type = GENE_NAME)

    # save data
    old_ANNOTATION_TYPE <<- ANNOTATION_TYPE
    ANNOTATION_TYPE <<- paste0(old_ANNOTATION_TYPE, '_from_', TRIMMING_PROB_TYPE, '_MHprob', LIGATION_MH_PARAM, '_trimMHprob', INT_MH_PARAM)
    if (NP_COND == TRUE){
        ANNOTATION_TYPE <<- paste0(ANNOTATION_TYPE, '_NPcond')
    }

    filename = processed_data_path(sample_annotation = SAMPLE_ANNOT)
    fwrite(processed, filename, sep = '\t')
}
