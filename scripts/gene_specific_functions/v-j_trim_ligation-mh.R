stopifnot(LOCUS == 'alpha')
stopifnot(INSERTIONS == 'zero')

get_oriented_whole_nucseqs <- function(whole_nucseq = get_whole_nucseqs()){
    # reorient sequence so that it is 5 -> 3 on the actual trimmed strand
    whole_nucseq[substring(gene, 4, 4) == 'J', sequence := unlist(lapply(sequence, function(x) as.character(reverseComplement(DNAString(x)))))]
    whole_nucseq[substring(gene, 4, 4) == 'J', j_gene_sequence := sequence]
    whole_nucseq[substring(gene, 4, 4) == 'V', v_gene_sequence := sequence]
    return(whole_nucseq[, -c('sequence')])
}

source(paste0(MOD_PROJECT_PATH,'/scripts/gene_count_specific_functions/double.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/mh_functions.R'))

filter_motif_data_for_possible_sites <- function(motif_data, whole_nucseq = get_oriented_whole_nucseqs(), gene_type = GENE_NAME, trim_type = TRIM_TYPE, only_nonprod_sites = ONLY_NONPROD_SITES){
    # Ensure that insertions are set to 'zero'
    stopifnot(INSERTIONS == 'zero')

    # Define frame-related columns
    frame_cols = c('frame_type', 'frame_stop')

    # Remove frame-related columns if present in motif_data
    if (any(frame_cols %in% colnames(motif_data))){
        cols = colnames(motif_data)[!(colnames(motif_data) %in% frame_cols)]
        motif_data = motif_data[, ..cols]
    }

    # Retrieve gene order based on gene_type
    genes = get_gene_order(gene_type)

    # Read frame data and filter for possible sites
    possible_sites_subset = get_all_possible_sites(gene_type)

    # Merge motif data with possible sites subset
    cols2 = c(paste0(genes), 'v_trim', 'j_trim', 'ligation_mh', 'processed_sequence')
    cols3 = unique(c(cols2, colnames(motif_data)[colnames(motif_data) %in% colnames(possible_sites_subset)]))
    tog = merge(motif_data, possible_sites_subset, by = cols3)

    # fill in remaining unobserved, but possible sites 
    ## Note: we are not including annotations that can be adjusted to a ligation-mh scenario; while these sites are possible, we are ignoring them since we have moved their counts to the ligation-mh scenario
    
    filled_tog = fill_in_missing_possible_sites(possible_sites_subset, tog, trim_type, gene_type)
    filled_tog[, index := .GRP, by = .(v_gene, j_gene, processed_sequence)]

    # remove counts for productive seqs, if noted
    if (only_nonprod_sites == TRUE){
        filled_tog[frame_type == 'In' & frame_stop == FALSE, count := NA]
    }
    return(filled_tog)
}

get_all_possible_sites <- function(gene_type = GENE_NAME){
    # Retrieve gene order based on gene_type
    genes = get_gene_order(gene_type)

    # Read frame data and filter for possible sites
    frame_data = read_frames_data()

    possible_sites = possible_sites[v_trim <= UPPER_TRIM_BOUND & j_trim <= UPPER_TRIM_BOUND]

    # Define columns for filtering possible sites
    cols = c(paste0(genes), 'frame_type', 'frame_stop', 'v_trim', 'j_trim', 'ligation_mh')
    cols = c(cols, 'processed_sequence', 'processed_nt_change', 'ligation_mh_nt')
    possible_sites_subset = unique(possible_sites[, ..cols])
    return(possible_sites_subset)
}

get_missing_possible_sites <- function(possible_sites, filtered_motif_data, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    genes = get_gene_order(gene_type)
    trims = get_trim_order(trim_type)

    # Get observed sets of gene pairs
    possible_sites$gene_pair = paste0(possible_sites$v_gene, '_', possible_sites$j_gene)
    filtered_motif_data$gene_pair = paste0(filtered_motif_data$v_gene, '_', filtered_motif_data$j_gene) 
    possible_sites_subset = possible_sites[gene_pair %in% unique(filtered_motif_data$gene_pair)]

    # get unobserved scenarios
    # cols = colnames(possible_sites_subset)[colnames(possible_sites_subset) %in% colnames(filtered_motif_data)]
    # unobserved = possible_sites_subset[!filtered_motif_data, on = cols]
    unobserved = possible_sites_subset[!filtered_motif_data, on = colnames(possible_sites_subset)]

    return(unobserved)
}

fill_in_missing_possible_sites <- function(possible_sites, filtered_motif_data, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    trims = get_trim_order(trim_type)
    genes = get_gene_order(gene_type)

    # get unobserved scenarios
    unobserved = get_missing_possible_sites(possible_sites, filtered_motif_data, trim_type, gene_type)

    if (nrow(unobserved) == 0){
        return(filtered_motif_data)
    } else {
        unobserved = get_oriented_full_sequences(unobserved, gene_type = gene_type)
        setnames(unobserved, paste0(genes, '_sequence'), paste0(genes, '_whole_seq'))

        unobserved = apply_get_nuc_context(unobserved, trim_type)

        unobserved$subject = unique(filtered_motif_data$subject)
        unobserved$gene_type = unique(filtered_motif_data$gene_type)
        unobserved[[paste0(trim_type, '_observed')]] = FALSE
        unobserved$count = 0

        together = rbind(filtered_motif_data, unobserved, fill = TRUE)
        return(together)
    }
}

get_igor_trimming_probs <- function(){
    jtrim_path = paste0(MOD_OUTPUT_PATH, '/meta_data/igor_alpha_j_trim_params.tsv')
    vtrim_path = paste0(MOD_OUTPUT_PATH, '/meta_data/igor_alpha_v_trim_params.tsv')

    if (!file.exists(vtrim_path)){
        print(paste0('need to produce trimming prob file'))
        py_script = paste0(PROJECT_PATH, '/mh_simulation_scripts/ligation-mh_signal_simulator/get_igor_params.py')
        command = paste('python', py_script, MOD_OUTPUT_PATH)
        
        system(command)
    }

    jt = fread(jtrim_path)
    vt = fread(vtrim_path)

    # renormalize
    vt[v_trim_prob == 0, v_trim_prob := 1e-30]
    jt[j_trim_prob == 0, j_trim_prob := 1e-30]

    # add additional trims
    vtrims = data.table(expand.grid(v_gene = unique(vt$v_gene),
                                    v_trim = seq(max(vt$v_trim) + 1, 20),
                                    v_trim_prob = 1e-30))
    jtrims = data.table(expand.grid(j_gene = unique(jt$j_gene),
                                    j_trim = seq(max(jt$j_trim) + 1, 20),
                                    j_trim_prob = 1e-30))

    vt = rbind(vt, vtrims)
    jt = rbind(jt, jtrims)

    vt[, v_trim_prob := v_trim_prob/sum(v_trim_prob), by = v_gene]
    jt[, j_trim_prob := j_trim_prob/sum(j_trim_prob), by = j_gene]

    return(list(v_trim = vt, j_trim = jt))
}

get_all_annotations <- function(tcr_dataframe, configs = read_frames_data()){
    cols = c('v_gene', 'j_gene', 'v_trim', 'j_trim', 'processed_sequence', 'ligation_mh')
    if ('productive' %in% colnames(tcr_dataframe)){
        cols = c(cols, 'productive')
    }
    if ('vj_insert' %in% colnames(tcr_dataframe)){
        cols = c(cols, 'vj_insert')
    }
    if ('count' %in% colnames(tcr_dataframe)){
        cols = c(cols, 'count')
    } else {
        tcr_dataframe$count = 1
        cols = c(cols, 'count')
    }
    if (!('ligation_mh' %in% colnames(tcr_dataframe))){
        tcr_dataframe$ligation_mh = 0
    }
    tcr_dataframe_wide = merge(tcr_dataframe, configs, by = c('v_gene', 'j_gene', 'v_trim', 'j_trim', 'ligation_mh'))
    tcr_dataframe_wide = tcr_dataframe_wide[, ..cols]
    tcr_dataframe_wide[, index := 1:.N]
    setnames(tcr_dataframe_wide, 'v_trim', 'original_v_trim')
    setnames(tcr_dataframe_wide, 'j_trim', 'original_j_trim')
    setnames(tcr_dataframe_wide, 'ligation_mh', 'original_ligation_mh')

    # combine by processed sequence
    cols = c('v_gene', 'j_gene', 'processed_sequence', 'vj_insert')
    condensed = tcr_dataframe_wide[, sum(count), by = cols]
    setnames(condensed, 'V1', 'count')
    condensed[, index := 1:.N]

    # get all other possible annotations
    all_annotations = merge(configs, condensed, by = c('v_gene', 'j_gene', 'processed_sequence'), allow.cartesian = TRUE)

    return(all_annotations)
}

adjust_trimming_sites_for_ligation_mh <- function(tcr_dataframe, sample_annotation){
    configs = read_frames_data()
    trimming_probs = get_igor_trimming_probs()

    cols = c('v_gene', 'j_gene', 'processed_sequence', 'v_trim', 'j_trim', 'ligation_mh', 'vj_insert', 'count')
    configs = merge(configs, trimming_probs$v_trim[, -c('prob_sum')], by = c('v_gene', 'v_trim'))
    configs = merge(configs, trimming_probs$j_trim[, -c('prob_sum')], by = c('j_gene', 'j_trim'))

    configs[, joint_trimming_prob := v_trim_prob * j_trim_prob]
    configs[joint_trimming_prob == 0, joint_trimming_prob := 1e-20]

    # get all other possible annotations
    all_annotations = get_all_annotations(tcr_dataframe, configs)

    if (sample_annotation){
        # sample annotations
        set.seed(123)
        sampled_annotations = all_annotations[, .SD[sample(.N, 1, prob = joint_trimming_prob)], by = index]

        TRIMMING_LIGATION_REANNOTATED <<- TRUE
        cols = c('index', cols)
    } else {
        cols = c('index', cols)
        sampled_annotations = all_annotations
    }
    return(sampled_annotations[, ..cols])
}
