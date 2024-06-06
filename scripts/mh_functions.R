get_pnucs <- function(whole_gene_nucseq, orient, pnuc_count){
    whole_gene_nucseq = DNAStringSet(whole_gene_nucseq)
    stopifnot(orient %in% c('top', 'bottom'))
    # top means top strand oriented 5 > 3
    # bottom means bottom strand oriented 3 > 5
    if (orient == 'top'){
        possible_pnucs = substring(as.character(reverseComplement(whole_gene_nucseq)),1, pnuc_count) 
    } else {
        possible_pnucs = substring(as.character(reverseComplement(whole_gene_nucseq)), nchar(whole_gene_nucseq) - pnuc_count + 1, nchar(whole_gene_nucseq))
    }
    return(possible_pnucs)
}

get_overlapping_regions_ligation_mh <- function(v_gene_top_seq, j_gene_bottom_seq, v_trim, j_trim, pnucs = 2){
    require(Biostrings)
    v_pnucs = get_pnucs(v_gene_top_seq, 'top', pnucs)    
    j_pnucs = get_pnucs(j_gene_bottom_seq, 'bottom', pnucs)    
    
    v_gene_top_seq_p = paste0(v_gene_top_seq, v_pnucs)
    j_gene_bottom_seq_p = paste0(j_pnucs, j_gene_bottom_seq)
    
    # get overlapping seqs when aligning trim sites
    v_overlap = substring(v_gene_top_seq_p, nchar(v_gene_top_seq_p) - (2*pnucs + v_trim + j_trim) + 1 , nchar(v_gene_top_seq_p))
    j_overlap = substring(j_gene_bottom_seq_p, 1, (2*pnucs + v_trim + j_trim))

    # get v_gene.j_trimmed overlaps
    vg_jt = substring(v_overlap, 1, pnucs+j_trim)
    jt_vg = substring(j_overlap, 1, pnucs+j_trim)

    # get j_gene.v_trimmed overlaps
    jg_vt = substring(j_overlap, pnucs+j_trim + 1)
    vt_jg = substring(v_overlap, pnucs+j_trim + 1)

    return(data.table('v_gene.j_trimmed' = vg_jt, 'j_trimmed.v_gene' = jt_vg, 'j_gene.v_trimmed' = jg_vt, 'v_trimmed.j_gene' = vt_jg))
}

get_mh_ligation_mh <- function(seq1, seq2, aligning_trim){
    stopifnot(aligning_trim %in% c('j_trim', 'v_trim'))
    stopifnot(all(nchar(seq1) == nchar(seq2)))
    
    compl = c('A' = 'T', 'T' = 'A', 'G' = 'C', 'C' = 'G')
    match = c('A' = 'A/T', 'T' = 'A/T', 'C' = 'G/C', 'G' = 'G/C')
    
    max_len = nchar(seq1)
    
    if (max_len == 0 | nchar(seq2) == 0){
        mh = data.table()
    } else {
        names = get_mh_colnames(aligning_trim, max_len)
        
        # Precompute constant values
        seq1_compl = compl[as.character(strsplit(seq1, NULL)[[1]])]
        seq1_match = match[as.character(strsplit(seq1, NULL)[[1]])]
        
        # Use vectorized operations to create the matrix
        mh_matrix = matrix('-', ncol = max_len)
        mh_matrix[, seq1_compl == strsplit(seq2, NULL)[[1]]] = seq1_match[seq1_compl == strsplit(seq2, NULL)[[1]]]
        
        # Convert the matrix to data.table
        mh = as.data.table(mh_matrix)
        colnames(mh) = names
    }
    return(mh)
}

get_mh_colnames <- function(aligning_trim, positions){
    if (aligning_trim == 'j_trim'){
        names = paste0(aligning_trim, '_mh_pos_', seq(positions, 1)) 
    } else {
        names = paste0(aligning_trim, '_mh_pos_', seq(1, positions)) 
    }
    return(names)
}

fill_in_missing_mh_positions <- function(mh_dt, max_len, aligning_trim){
    stopifnot(max_len >= ncol(mh_dt))

    all_names = get_mh_colnames(aligning_trim, max_len)
    not_present = all_names[!(all_names %in% colnames(mh_dt))]
    if (length(not_present) > 0){
        additional = as.data.table(matrix(ncol = length(not_present))) 
        colnames(additional) = not_present
        mh_dt = cbind(mh_dt, additional)
    }
    setcolorder(mh_dt, all_names)
    return(mh_dt)
}

get_mh_and_fill <- function(seq1, seq2, aligning_trim, max_len){
    return(fill_in_missing_mh_positions(get_mh_ligation_mh(seq1, seq2, aligning_trim), max_len, aligning_trim)) 
}

get_mh_dataframe <- function(data, aligning_trim, aligning_gene){
    col1 = paste0(aligning_trim, 'med.', aligning_gene)
    col2 = paste0(aligning_gene, '.', aligning_trim, 'med')
    stopifnot(col1 %in% colnames(data))
    stopifnot(col2 %in% colnames(data))

    cols = c(col1, col2)
    subset = unique(data[, ..cols]) 
    max = max(nchar(data[[col1]]), nchar(data[[col2]]))
    # max = UPPER_TRIM_BOUND + PNUC_COUNT
    names = get_mh_colnames(aligning_trim, max)

    if (nrow(subset) == 1){
        subset[, paste0(names) := data.table(matrix(mapply(get_mh_and_fill, get(col1), get(col2), aligning_trim, max), nrow = 1))]
    } else {
        subset[, paste0(names) := as.data.table(apply(t(mapply(get_mh_and_fill, get(col1), get(col2), aligning_trim, max)), 2, unlist))]
    }

    data = merge(data, subset, by = cols)
    return(data)
}

get_mh <- function(seq1, seq2){
    seq2_comp = as.character(complement(DNAStringSet(seq2)))
    mh = mcmapply(function(X,Y) sum(str_count(X,Y)), strsplit(seq1, ''), strsplit(seq2_comp, ''))
    return(mh)
}

get_fully_contiguous_mh <- function(seq1, seq2){
    noncontig_mh = get_mh(seq1, seq2)
    contig_mh = ifelse(noncontig_mh == nchar(seq1), noncontig_mh, 0)
    return(contig_mh)
}

get_overlapping_regions <- function(v_gene_top_seq, j_gene_bottom_seq, v_trim, j_trim, overlap_count, pnucs = 2, positions = c('up', 'mid', 'down'), rename_cols = TRUE){
    # j_gene_bottom_seq is the bottom strand oriented 5 > 3
    # this reorientation will orient the bottom strand 3 > 5 so that it is complementary to the V-gene top strand
    j_gene_bottom_seq = reorient_j_bottom_strand(j_gene_bottom_seq)

    require(Biostrings)
    v_pnucs = get_pnucs(v_gene_top_seq, 'top', pnucs)    
    j_pnucs = get_pnucs(j_gene_bottom_seq, 'bottom', pnucs)    
    
    v_gene_top_seq_p = paste0(v_gene_top_seq, v_pnucs)
    j_gene_bottom_seq_p = paste0(j_pnucs, j_gene_bottom_seq)
    
    # get overlapping seqs when aligning trim sites
    v_overlap = substring(v_gene_top_seq_p, nchar(v_gene_top_seq_p) - (2*pnucs + v_trim + j_trim + overlap_count) + 1 , nchar(v_gene_top_seq_p))
    j_overlap = substring(j_gene_bottom_seq_p, 1, (2*pnucs + v_trim + j_trim + overlap_count))

    # get v_gene.j_trimmed overlaps
    vg_jt = substring(v_overlap, 1, pnucs+j_trim)
    jt_vg = substring(j_overlap, 1, pnucs+j_trim)

    # get j_gene.v_trimmed overlaps
    jg_vt = substring(j_overlap, pnucs+j_trim + 1 + overlap_count)
    vt_jg = substring(v_overlap, pnucs+j_trim + 1 + overlap_count)

    # get overlaps
    if (length(overlap_count) == 1){
        overlap_count_long = rep(overlap_count, length(v_overlap))
    } else {
        overlap_count_long = overlap_count
    }

    v_mid = ifelse(overlap_count_long > 0, substring(v_overlap, pnucs + j_trim + 1, pnucs + j_trim + overlap_count), '')
    j_mid = ifelse(overlap_count_long > 0, substring(j_overlap, pnucs + j_trim + 1, pnucs + j_trim + overlap_count), '')

    up = data.table(vg_jt, jt_vg)
    down = data.table(vt_jg, jg_vt)
    mid = data.table(v_mid, j_mid)

    overlaps = c()
    for (pos in positions){
        overlaps = cbind(overlaps, get(pos))
    }

    if (rename_cols){
        names = get_overlap_names(overlap_count, positions)
        setnames(overlaps, colnames(overlaps), names)
    }
    return(overlaps)
}

get_processed_sequences <- function(adjusted_df, frame_data, pnucs = 2){
    stopifnot('ligation_mh' %in% colnames(adjusted_df))
    # j_gene_sequence should be the bottom strand oriented 5>3
    cols = colnames(adjusted_df)
    # get overlapping sequences
    overlaps = adjusted_df[, get_overlapping_regions(v_gene_sequence, j_gene_sequence, v_trim, j_trim, ligation_mh, positions = c('up', 'mid', 'down'), rename_cols = FALSE)]
    tog = cbind(adjusted_df, overlaps) 

    # get bottom gene sequence 
    j_gene_bottom_seq = reorient_j_bottom_strand(adjusted_df$j_gene_sequence)

    # get pnucs
    require(Biostrings)
    v_pnucs = get_pnucs(adjusted_df$v_gene_sequence, 'top', pnucs)    
    j_pnucs = get_pnucs(j_gene_bottom_seq, 'bottom', pnucs)    

    v_gene_top_seq_p = paste0(adjusted_df$v_gene_sequence, v_pnucs)
    j_gene_bottom_seq_p = paste0(j_pnucs, j_gene_bottom_seq)

    j_gene_top_seq_p = as.character(complement(DNAStringSet(j_gene_bottom_seq_p)))

    tog$v_gene_sequence_pnucs = v_gene_top_seq_p
    tog$j_gene_sequence_pnucs_top = j_gene_top_seq_p

    # trim sequences
    tog[, trimmed_v_gene_sequence := substring(v_gene_sequence_pnucs, 1, nchar(v_gene_sequence_pnucs) - 2 - v_trim)]
    tog[, trimmed_j_gene_sequence := substring(j_gene_sequence_pnucs_top, j_trim + 3)]

    # check MH
    tog[ligation_mh > 0, v_gene_mh := substring(trimmed_v_gene_sequence, nchar(trimmed_v_gene_sequence) - ligation_mh + 1)]
    tog[ligation_mh > 0, j_gene_mh := substring(trimmed_j_gene_sequence, 1, ligation_mh)]
    stopifnot(all(tog[ligation_mh > 0]$v_mid == tog[ligation_mh > 0]$v_gene_mh))
    stopifnot(all(tog[ligation_mh > 0]$v_mid == tog[ligation_mh > 0]$j_gene_mh))

    # get sequences corrected for MH
    tog[, trimmed_j_gene_sequence_wo_mh := substring(trimmed_j_gene_sequence, ligation_mh + 1)]

    # get processed sequences
    tog[, processed_sequence := paste(trimmed_v_gene_sequence, trimmed_j_gene_sequence_wo_mh, sep = '')]
    tog[, processed_nt_change := v_trim + j_trim + ligation_mh]

    # get cdr3 sequences
    cindex = get_cdr3_indices(frame_data)
    tog = merge(tog, cindex$v, by = 'v_gene')
    tog = merge(tog, cindex$j, by = 'j_gene')

    tog[, processed_protein := as.character(translate(DNAStringSet(processed_sequence)))]

    tog[, processed_cdr3 := substring(processed_protein, cys_index, nchar(processed_protein) - phe_from_end_index)]

    setnames(tog, 'v_mid', 'ligation_mh_nt')
    cols = c(cols, 'processed_sequence', 'processed_nt_change', 'ligation_mh_nt', 'processed_cdr3', 'cys_character', 'phe_character')
    return(tog[, ..cols])
}

get_cdr3_indices <- function(frame_data){
    vframe = frame_data[id %like% 'TRAV']
    vframe[, c("extra", "cdr1", "cdr2", 'cdr3') := tstrsplit(cdr_columns, ";", fixed = TRUE)]
    vframe[, c("cdr3_start", "CDR3_end") := tstrsplit(cdr3, "-", fixed = TRUE)]
    vframe[, seq_start := substring(aligned_protseq, 1, cdr3_start)]
    vframe[, dot_count := stringr::str_count(seq_start, "\\.")]
    vframe[, cys_index := as.numeric(cdr3_start) - dot_count]
    vframe[, cys_character := substring(seq_start, nchar(seq_start))]
    vframe[cys_character == '.', cys_character := ""]
    setnames(vframe, 'id', 'v_gene')

    jframe = frame_data[id %like% 'TRAJ']
    jframe[, c("cdr3_start", "CDR3_end") := tstrsplit(cdr_columns, "-", fixed = TRUE)]
    jframe[, seq_start := substring(aligned_protseq, 1, CDR3_end)]
    jframe[, seq_end := substring(aligned_protseq, CDR3_end)]
    jframe[, dot_count := stringr::str_count(seq_start, "\\.")]
    jframe[, dot_count_end := stringr::str_count(seq_end, "\\.")]
    jframe[, phe_index := as.numeric(CDR3_end) - dot_count]
    jframe[, phe_from_end_index := nchar(seq_end) - 1 - dot_count_end]
    jframe[, phe_character := substring(seq_end, 1, 1)]
    setnames(jframe, 'id', 'j_gene')

    return(list(v = vframe[, c('v_gene', 'cys_index', 'cys_character')], j = jframe[, c('j_gene', 'phe_index', 'phe_from_end_index', 'phe_character')]))
}

get_possible_ligation_mh_fixed_trim <- function(data, overlap_count){
    # get overlapping regions
    names = get_overlap_names(overlap_count, 'mid')
    subset = data[, get_overlapping_regions(v_gene_sequence, j_gene_sequence, v_trim, j_trim, overlap_count, positions = 'mid')]
    colnames(subset) = names

    # get MH
    pos = 'mid'
    v_seq_col = paste0('v_gene_', pos, '_overlap_', overlap_count)
    j_seq_col = paste0('j_gene_', pos, '_overlap_', overlap_count)

    n = get_fully_contiguous_mh(subset[[v_seq_col]], subset[[j_seq_col]])
    return(n)
}
