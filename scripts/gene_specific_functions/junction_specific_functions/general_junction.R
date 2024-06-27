get_overlapping_regions <- function(top_gene_name, bottom_gene_name, top_gene_top_seq, bottom_gene_bottom_seq, top_trim, bottom_trim, overlap_count, pnucs = 2, positions = c('up', 'mid', 'down'), rename_cols = TRUE){
    # bottom_gene_bottom_seq is the bottom strand oriented 5 > 3
    # this reorientation will orient the bottom strand 3 > 5 so that it is
    # complementary to the top-gene top strand
    bottom_gene_bottom_seq = reorient_j_bottom_strand(bottom_gene_bottom_seq)

    require(Biostrings)
    top_pnucs = get_pnucs(top_gene_top_seq, 'top', pnucs)    
    bottom_pnucs = get_pnucs(bottom_gene_bottom_seq, 'bottom', pnucs)    
    
    top_gene_top_seq_p = paste0(top_gene_top_seq, top_pnucs)
    bottom_gene_bottom_seq_p = paste0(bottom_pnucs, bottom_gene_bottom_seq)
    
    # get overlapping seqs when aligning trim sites
    top_overlap = substring(top_gene_top_seq_p, nchar(top_gene_top_seq_p) - (2*pnucs + top_trim + bottom_trim + overlap_count) + 1 , nchar(top_gene_top_seq_p))
    bottom_overlap = substring(bottom_gene_bottom_seq_p, 1, (2*pnucs + top_trim + bottom_trim + overlap_count))

    # get top_gene.bottom_trimmed overlaps
    topg_bottomt = substring(top_overlap, 1, pnucs+bottom_trim)
    bottomt_topg = substring(bottom_overlap, 1, pnucs+bottom_trim)

    # get bottom_gene.top_trimmed overlaps
    bottomg_topt = substring(bottom_overlap, pnucs+bottom_trim + 1 + overlap_count)
    topt_bottomg = substring(top_overlap, pnucs+bottom_trim + 1 + overlap_count)

    # get overlaps
    if (length(overlap_count) == 1){
        overlap_count_long = rep(overlap_count, length(top_overlap))
    } else {
        overlap_count_long = overlap_count
    }

    top_mid = ifelse(overlap_count_long > 0, substring(top_overlap, pnucs + bottom_trim + 1, pnucs + bottom_trim + overlap_count), '')
    bottom_mid = ifelse(overlap_count_long > 0, substring(bottom_overlap, pnucs + bottom_trim + 1, pnucs + bottom_trim + overlap_count), '')

    up = data.table(topg_bottomt, bottomt_topg)
    down = data.table(topt_bottomg, bottomg_topt)
    mid = data.table(top_mid, bottom_mid)

    overlaps = c()
    for (pos in positions){
        overlaps = cbind(overlaps, get(pos))
    }

    if (rename_cols){
        names = get_overlap_names(overlap_count, positions, top_gene_name, bottom_gene_name)
        setnames(overlaps, colnames(overlaps), names)
    }
    return(overlaps)
}

get_overlapping_regions_ligation_mh <- function(top_gene_top_seq, bottom_gene_bottom_seq, top_trim, bottom_trim, pnucs = 2){
    require(Biostrings)
    top_pnucs = get_pnucs(top_gene_top_seq, 'top', pnucs)    
    bottom_pnucs = get_pnucs(bottom_gene_bottom_seq, 'bottom', pnucs)    
    
    top_gene_top_seq_p = paste0(top_gene_top_seq, top_pnucs)
    bottom_gene_bottom_seq_p = paste0(bottom_pnucs, bottom_gene_bottom_seq)
    
    # get overlapping seqs when aligning trim sites
    top_otoperlap = substring(top_gene_top_seq_p, nchar(top_gene_top_seq_p) - (2*pnucs + top_trim + bottom_trim) + 1 , nchar(top_gene_top_seq_p))
    bottom_overlap = substring(bottom_gene_bottom_seq_p, 1, (2*pnucs + top_trim + bottom_trim))

    # get top_gene.bottom_trimmed overlaps
    topg_bottomt = substring(top_overlap, 1, pnucs+bottom_trim)
    bottomt_topg = substring(bottom_overlap, 1, pnucs+bottom_trim)

    # get bottom_gene.top_trimmed overlaps
    bottomg_topt = substring(bottom_overlap, pnucs+bottom_trim + 1)
    topt_bottomg = substring(top_overlap, pnucs+bottom_trim + 1)

    return(data.table('top_gene.bottom_trimmed' = topg_bottomt, 'bottom_trimmed.top_gene' = bottomt_topg, 'bottom_gene.top_trimmed' = bottomg_topt, 'top_trimmed.bottom_gene' = topt_bottomg))
}

get_possible_ligation_mh_fixed_trim <- function(data, overlap_count, top_gene_sequence_name, bottom_gene_sequence_name, top_sequence_trim_name, bottom_sequence_trim_name){
    top_gene = substring(top_gene_sequence_name, 1, 6)
    bottom_gene = substring(bottom_gene_sequence_name, 1, 6)

    # get overlapping regions
    names = get_overlap_names(overlap_count, 'mid', top_gene, bottom_gene)
    subset = data[, get_overlapping_regions(top_gene, bottom_gene, get(top_gene_sequence_name), get(bottom_gene_sequence_name), get(top_sequence_trim_name), get(bottom_sequence_trim_name), overlap_count, positions = 'mid')]
    colnames(subset) = names

    # get MH
    pos = 'mid'
    v_seq_col = paste0(top_gene, '_', pos, '_overlap_', overlap_count)
    j_seq_col = paste0(bottom_gene, '_', pos, '_overlap_', overlap_count)

    n = get_fully_contiguous_mh(subset[[v_seq_col]], subset[[j_seq_col]])
    return(n)
}

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
    trims = get_trim_vars(trim_type)

    # Read frame data and filter for possible sites
    possible_sites_subset = get_all_possible_sites(gene_type)

    # Merge motif data with possible sites subset
    cols2 = c(genes, trims, 'processed_sequence')
    cols3 = unique(c(cols2, colnames(motif_data)[colnames(motif_data) %in% colnames(possible_sites_subset)]))
    tog = merge(motif_data, possible_sites_subset, by = cols3)

    # fill in remaining unobserved, but possible sites 
    ## Note: we are not including annotations that can be adjusted to a ligation-mh scenario; while these sites are possible, we are ignoring them since we have moved their counts to the ligation-mh scenario
    
    filled_tog = fill_in_missing_possible_sites(possible_sites_subset, tog, trim_type, gene_type)
    cols = c(genes, 'processed_sequence')
    filled_tog[, index := .GRP, by = cols]

    # remove counts for productive seqs, if noted
    if (only_nonprod_sites == TRUE){
        if (grepl('nonprod', PARAM_GROUP)){
            filled_tog[frame_type == 'In' & frame_stop == FALSE, count := NA]
        } else if (grepl('prod', PARAM_GROUP)){
            filled_tog[frame_type == 'Out' | frame_stop == TRUE, count := NA]
        }
    }
    return(filled_tog)
}

get_all_possible_sites <- function(gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    # Retrieve gene order based on gene_type
    genes = get_gene_order(gene_type)
    trims = get_trim_vars(trim_type)

    # Read frame data and filter for possible sites
    possible_sites = read_frames_data()

    possible_sites = possible_sites[get(trims[1]) <= UPPER_TRIM_BOUND]
    possible_sites = possible_sites[get(trims[2]) <= UPPER_TRIM_BOUND]

    # Define columns for filtering possible sites
    cols = c(paste0(genes), 'frame_type', 'frame_stop', trims, 'processed_sequence', 'processed_nt_change')
    possible_sites_subset = unique(possible_sites[, ..cols])
    return(possible_sites_subset)
}

get_missing_possible_sites <- function(possible_sites, filtered_motif_data, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    genes = get_gene_order(gene_type)
    trims = get_trim_order(trim_type)

    # Get observed sets of gene pairs
    possible_sites$gene_pair = paste0(possible_sites[[genes[1]]], '_', possible_sites[[genes[2]]])
    filtered_motif_data$gene_pair = paste0(filtered_motif_data[[genes[1]]], '_', filtered_motif_data[[genes[2]]]) 
    possible_sites_subset = possible_sites[gene_pair %in% unique(filtered_motif_data$gene_pair)]

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

adjust_trimming_sites_for_ligation_mh <- function(tcr_dataframe, sample_annotation, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    stopifnot(sample_annotation == FALSE)
    configs = read_frames_data()
    # get all other possible annotations
    all_annotations = get_all_annotations(tcr_dataframe, configs, trim_type, gene_type)
    return(all_annotations)
}

get_all_annotations <- function(tcr_dataframe, configs = read_frames_data(), trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    genes = get_gene_order(gene_type)
    trims = get_trim_vars(trim_type)

    cols = c(genes, trims, 'processed_sequence')
    if ('productive' %in% colnames(tcr_dataframe)){
        cols = c(cols, 'productive')
    }
    if (JOINING_INSERT %in% colnames(tcr_dataframe)){
        cols = c(cols, JOINING_INSERT)
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
    cols2 = c(genes, trims)
    tcr_dataframe_wide = merge(tcr_dataframe, configs, by = cols2)
    tcr_dataframe_wide = tcr_dataframe_wide[, ..cols]
    tcr_dataframe_wide[, index := 1:.N]
    setnames(tcr_dataframe_wide, trims[1], paste0('original_', trims[1]))
    setnames(tcr_dataframe_wide, trims[2], paste0('original_', trims[2]))
    setnames(tcr_dataframe_wide, 'ligation_mh', 'original_ligation_mh')

    # combine by processed sequence
    cols = c(genes, 'processed_sequence', JOINING_INSERT)
    condensed = tcr_dataframe_wide[, sum(count), by = cols]
    setnames(condensed, 'V1', 'count')
    condensed[, index := 1:.N]

    # get all other possible annotations
    all_annotations = merge(configs, condensed, by = c(genes, 'processed_sequence'), allow.cartesian = TRUE)

    return(all_annotations)
}
