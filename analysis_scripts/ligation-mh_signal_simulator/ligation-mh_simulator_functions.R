get_igor_gene_usage_params <- function(){
    # get igor params
    j_choice_path= paste0(MOD_OUTPUT_PATH, '/meta_data/igor_alpha_jchoice_params.tsv')
    v_choice_path= paste0(MOD_OUTPUT_PATH, '/meta_data/igor_alpha_vchoice_params.tsv')

    jg = fread(j_choice_path)
    vg = fread(v_choice_path)
    return(list(v_choice = vg, j_choice = jg))
}

get_trimming_probs <- function(type, int_mh_prob, configs){
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

    # adjust trimming probabilities by MH param
    trimming_probs = adjust_trimming_probs_by_mh(vt, jt, configs, int_mh_prob)
    return(trimming_probs)
}

get_ligation_probabilities <- function(lig_param, possible_ligs){
    stopifnot(lig_param >= 0)
    intercept = 1
    lig_probs = possible_ligs*lig_param + intercept
    names(lig_probs) = possible_ligs    
    lig_probs['no_ligation'] = 0.1
    prob_sum = sum(lig_probs)
    lig_probs = lig_probs/prob_sum
    return(lig_probs)
}

adjust_trimming_probs_by_mh <- function(v_probs, j_probs, configs, int_mh_prob){
    stopifnot(int_mh_prob >= 0)

    configs = get_mh_config_count(configs)
    configs = merge(configs, v_probs, by = c('v_gene', 'v_trim'))
    configs = merge(configs, j_probs, by = c('j_gene', 'j_trim'))

    cols = c('v_gene', 'j_gene', 'v_trim', 'j_trim', 'v_trim_prob', 'j_trim_prob', 'mh_config_count')

    combo_mh = unique(configs[, ..cols])

    combo_mh[, joint_trim_prob := v_trim_prob * j_trim_prob]

    combo_mh[, int_mh_covariate_function := int_mh_prob*mh_config_count]
    combo_mh[, int_mh_softmax_per_pair := exp(int_mh_covariate_function)/sum(exp(int_mh_covariate_function)), by = .(v_gene, j_gene)]

    combo_mh[, mh_adj_joint_trim_prob := joint_trim_prob * int_mh_softmax_per_pair]
    combo_mh[, mh_adj_joint_trim_prob := mh_adj_joint_trim_prob/sum(mh_adj_joint_trim_prob), by = .(v_gene, j_gene)]

    cols = c(cols, 'mh_adj_joint_trim_prob')

    combo_mh_subset = combo_mh[, ..cols]

    return(combo_mh_subset)
}
