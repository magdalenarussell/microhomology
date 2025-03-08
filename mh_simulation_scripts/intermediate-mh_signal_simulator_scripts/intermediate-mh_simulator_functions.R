get_igor_gene_usage_params <- function(){
    # get igor params
    j_choice_path= paste0(MOD_OUTPUT_PATH, '/meta_data/igor_alpha_jchoice_params.tsv')
    v_choice_path= paste0(MOD_OUTPUT_PATH, '/meta_data/igor_alpha_vchoice_params.tsv')

    jg = fread(j_choice_path)
    vg = fread(v_choice_path)
    return(list(v_choice = vg, j_choice = jg))
}

get_trimming_probs <- function(type){
    stopifnot(type %in% c('igor', 'motif_two-side-base-count-beyond', 'uniform', 'mh_adjusted_motif-two-side-base-count-beyond'))
    if (type == 'igor'){
        jtrim_path = paste0(MOD_OUTPUT_PATH, '/meta_data/igor_alpha_j_trim_params.tsv')
        vtrim_path = paste0(MOD_OUTPUT_PATH, '/meta_data/igor_alpha_v_trim_params.tsv')

        jt = fread(jtrim_path)
        vt = fread(vtrim_path)
    } else if (type == 'uniform'){
        jtrim_path = paste0(MOD_OUTPUT_PATH, '/meta_data/igor_alpha_j_trim_params.tsv')
        vtrim_path = paste0(MOD_OUTPUT_PATH, '/meta_data/igor_alpha_v_trim_params.tsv')

        jt = fread(jtrim_path)
        vt = fread(vtrim_path)

        vt[, v_trim_prob := 1]
        jt[, j_trim_prob := 1]
    } else if (type == 'motif_two-side-base-count-beyond'){
        # get predictions (i.e. probabilities) from no-MH model
        v_probs_path = get_model_predictions_file_path('False', model_type='motif_two-side-base-count-beyond')
        j_probs_path = get_model_predictions_file_path('False', model_type='motif_two-side-base-count-beyond')
        v_probs_path = str_replace(v_probs_path, ANNOTATION_TYPE, 'igor_sim_alpha')
        j_probs_path = str_replace(j_probs_path, ANNOTATION_TYPE, 'igor_sim_alpha')
        v_probs_path = str_replace(v_probs_path, PARAM_GROUP, 'nonproductive_v_trim')
        j_probs_path = str_replace(j_probs_path, PARAM_GROUP, 'nonproductive_j_trim')

        if (!file.exists(v_probs_path)){
            print(paste0('need to produce predictions file: ', v_probs_path, ' before trimming probabilities can be obtained'))
        }

        if (!file.exists(j_probs_path)){
            print(paste0('need to produce predictions file: ', j_probs_path, ' before trimming probabilities can be obtained'))
        }

        vt = fread(v_probs_path)
        vcols = c('v_gene', 'v_trim', 'predicted_prob')
        vt = unique(vt[, ..vcols])
        setnames(vt, 'predicted_prob', 'v_trim_prob')

        jt = fread(j_probs_path)
        jcols = c('j_gene', 'j_trim', 'predicted_prob')
        jt = unique(jt[, ..jcols])
        setnames(jt, 'predicted_prob', 'j_trim_prob')
    } else{
        # get predictions (i.e. probabilities) from no-MH model
        v_probs_path = get_model_predictions_file_path('False', model_type='motif_two-side-base-count-beyond')
        j_probs_path = get_model_predictions_file_path('False', model_type='motif_two-side-base-count-beyond')
        v_probs_path = str_replace(v_probs_path, ANNOTATION_TYPE, 'igor_sim_alpha')
        j_probs_path = str_replace(j_probs_path, ANNOTATION_TYPE, 'igor_sim_alpha')
        v_probs_path = str_replace(v_probs_path, PARAM_GROUP, 'nonproductive_v_trim_adjusted_mh')
        j_probs_path = str_replace(j_probs_path, PARAM_GROUP, 'nonproductive_j_trim_adjusted_mh')

        if (!file.exists(v_probs_path)){
            print(paste0('need to produce predictions file: ', v_probs_path, ' before trimming probabilities can be obtained'))
        }

        if (!file.exists(j_probs_path)){
            print(paste0('need to produce predictions file: ', j_probs_path, ' before trimming probabilities can be obtained'))
        }

        vt = fread(v_probs_path)
        vcols = c('v_gene', 'v_trim_adjusted_mh', 'predicted_prob')
        vt = unique(vt[, ..vcols])
        setnames(vt, 'predicted_prob', 'v_trim_prob')
        setnames(vt, 'v_trim_adjusted_mh', 'v_trim')

        jt = fread(j_probs_path)
        jcols = c('j_gene', 'j_trim_adjusted_mh', 'predicted_prob')
        jt = unique(jt[, ..jcols])
        setnames(jt, 'predicted_prob', 'j_trim_prob')
        setnames(jt, 'j_trim_adjusted_mh', 'j_trim')
    }

    vt = vt[v_trim >= LOWER_TRIM_BOUND & v_trim <= UPPER_TRIM_BOUND]
    jt = jt[j_trim >= LOWER_TRIM_BOUND & j_trim <= UPPER_TRIM_BOUND]
    return(list(v_trim = vt, j_trim = jt))
}

get_all_possible_configs <- function(gene_probs, trim_probs, int_mh_prob){
    stopifnot(int_mh_prob >= 0)

    source(paste0(MOD_PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/mh.R'))

    # get all gene and trimming combos
    jcols = c('j_gene', 'j_gene_sequence')
    vcols = c('v_gene', 'v_gene_sequence')
    
    jgene_seq = gene_probs$j_choice[, ..jcols]
    vgene_seq = gene_probs$v_choice[, ..vcols]

    jgene_seq_list = split(jgene_seq, seq(nrow(jgene_seq)))
    vgene_seq_list = split(vgene_seq, seq(nrow(vgene_seq)))

    trims = unique(trim_probs$v_trim$v_trim)

    combos = CJ(vgene_seq_list, jgene_seq_list, v_trim = trims, j_trim = trims, sorted = FALSE)

    combos[, c('v_gene', 'v_gene_sequence') := .(unlist(sapply(vgene_seq_list, function(dt) dt$v_gene)),
                                                 unlist(sapply(vgene_seq_list, function(dt) dt$v_gene_sequence)))]
    combos[, c('j_gene', 'j_gene_sequence') := .(unlist(sapply(jgene_seq_list, function(dt) dt$j_gene)),
                                                 unlist(sapply(jgene_seq_list, function(dt) dt$j_gene_sequence)))]

    cols = colnames(combos)[!(colnames(combos) %like% 'list')]
    combos = combos[, ..cols]

    combo_mh = process_for_mh(combos, whole_nucseq = get_oriented_whole_nucseqs(), overlap_vector = c(1, 2, 3, 4), trim_type = TRIM_TYPE, gene_type = GENE_NAME, prop = FALSE, positions = c('mid'))

    cols = c('v_gene', 'j_gene', 'v_trim', 'j_trim', 'mh_count_mid_overlap_1', 'mh_count_mid_overlap_2', 'mh_count_mid_overlap_3', 'mh_count_mid_overlap_4')
    combo_mh = combo_mh[, ..cols]

    combo_mh = merge(combo_mh, trim_probs$v_trim, by = c('v_gene', 'v_trim'))
    combo_mh = merge(combo_mh, trim_probs$j_trim, by = c('j_gene', 'j_trim'))
    combo_mh = merge(combo_mh, gene_probs$v_choice[, c('v_gene', 'v_gene_sequence')], by = 'v_gene')
    combo_mh = merge(combo_mh, gene_probs$j_choice[, c('j_gene', 'j_gene_sequence')], by = 'j_gene')

    combo_mh[, joint_trim_prob := v_trim_prob * j_trim_prob]

    combo_mh[, average_interior_mh := (mh_count_mid_overlap_1 + mh_count_mid_overlap_2 + mh_count_mid_overlap_3 + mh_count_mid_overlap_4)/4]

    combo_mh[, int_mh_covariate_function := int_mh_prob*average_interior_mh]
    combo_mh[, int_mh_softmax_per_pair := exp(int_mh_covariate_function)/sum(exp(int_mh_covariate_function)), by = .(v_gene, j_gene)]

    combo_mh[, mh_adj_joint_trim_prob := joint_trim_prob * int_mh_softmax_per_pair]
    combo_mh[, mh_adj_joint_trim_prob := mh_adj_joint_trim_prob/sum(mh_adj_joint_trim_prob), by = .(v_gene, j_gene)]

    cols = c(cols, 'v_gene_sequence', 'j_gene_sequence', 'mh_adj_joint_trim_prob')

    combo_mh_subset = combo_mh[, ..cols]

    return(combo_mh_subset)
}
