get_oriented_whole_nucseqs <- function(whole_nucseq = get_whole_nucseqs()){
    # reorient sequence so that it is 5 -> 3 on the actual trimmed strand
    whole_nucseq[substring(gene, 4, 4) == 'J', sequence := unlist(lapply(sequence, function(x) as.character(reverseComplement(DNAString(x)))))]
    whole_nucseq[substring(gene, 4, 4) == 'J', j_gene_sequence := sequence]
    whole_nucseq[substring(gene, 4, 4) == 'V', v_gene_sequence := sequence]
    return(whole_nucseq[, -c('sequence')])
}

source(paste0(MOD_PROJECT_PATH,'/scripts/gene_count_specific_functions/double.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/mh_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/gene_specific_functions/junction_specific_functions/general_junction.R'))

get_processed_sequences <- function(adjusted_df, frame_data, pnucs = 2){
    stopifnot('ligation_mh' %in% colnames(adjusted_df))
    # j_gene_sequence should be the bottom strand oriented 5>3
    cols = colnames(adjusted_df)
    # get overlapping sequences
    overlaps = adjusted_df[, get_overlapping_regions('v_gene', 'j_gene', v_gene_sequence, j_gene_sequence, v_trim, j_trim, ligation_mh, positions = c('up', 'mid', 'down'), rename_cols = FALSE)]
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
    stopifnot(all(tog[ligation_mh > 0]$top_mid == tog[ligation_mh > 0]$v_gene_mh))
    stopifnot(all(tog[ligation_mh > 0]$top_mid == tog[ligation_mh > 0]$j_gene_mh))

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

    setnames(tog, 'top_mid', 'ligation_mh_nt')
    cols = c(cols, 'processed_sequence', 'processed_nt_change', 'ligation_mh_nt', 'processed_cdr3', 'cys_character', 'phe_character')
    return(tog[, ..cols])
}

get_frames_data <- function(){
    # This function will only return frame data for annotations that are observable, after re-annotation to accommodate cases of maximal ligation-mh
    # Read only necessary columns and apply filter
    frames = fread('https://raw.githubusercontent.com/phbradley/conga/master/conga/tcrdist/db/combo_xcr.tsv')[organism %like% 'human' & chain == CHAIN_SUBTYPE][substring(id, 3, 3) == substring(CHAIN_TYPE, 3, 3)]

    # Combine operations to reduce redundancy
    v = frames[region == 'V', .(v_gene = id, v_frame = frame, v_seq = nucseq)]
    j = frames[region == 'J', .(j_gene = id, j_frame = frame, j_seq = nucseq)]

    # Merge operations
    v$dummy = 1
    j$dummy = 1
    gene_pairs = merge(v, j, by = 'dummy', allow.cartesian = TRUE)

    # Get all trimming sites
    trims = data.table(expand.grid(v_trim = seq(LOWER_TRIM_BOUND, 16), 
                                   j_trim = seq(LOWER_TRIM_BOUND, 18)))
    trims$dummy = 1

    # Merge genes and trims
    all = merge(gene_pairs, trims, by = 'dummy', allow.cartesian = TRUE)[, -c('dummy')]
    # Get sequence lengths and subset data to necessary columns
    all[, c('v_seq_len', 'j_seq_len') := .(nchar(v_seq), nchar(j_seq))]

    cols = c('v_gene', 'j_gene', 'v_frame', 'j_frame', 'v_seq_len', 'j_seq_len', 'v_trim', 'j_trim', 'v_gene_sequence', 'j_gene_sequence')
    all = get_oriented_full_sequences(all)
    all = unique(all[, ..cols])

    # get possible ligation mh
    lig_mat = matrix(0, nrow = nrow(all), ncol = length(seq(1, 15)))

    for (overlap in seq(ncol(lig_mat))){
        lig = get_possible_ligation_mh_fixed_trim(all, overlap_count = overlap, 'v_gene_sequence', 'j_gene_sequence', 'v_trim', 'j_trim')
        lig_mat[, overlap] = lig
        print(paste0('finished getting ligation configurations for ', overlap, ' MH overlap'))
    }

    # get unique vals 
    unique_mh_list = apply(lig_mat, 1, unique)

    # combine
    adjusted_all = data.table(all, ligation_mh = unique_mh_list)

    # Expand the rows
    adjusted_grouped = as.data.table(adjusted_all[, unnest(.SD, cols = c("ligation_mh"))])

    # Get oriented full sequences and group genes by common features, also subset columns again
    cols2 = c(cols, 'ligation_mh')
    adjusted_grouped = unique(adjusted_grouped[, ..cols2])
    
    # get processed sequences
    adjusted_grouped = get_processed_sequences(adjusted_grouped, frames)

    # Get stop positions and frames
    adjusted_grouped = get_stop_codons(adjusted_grouped, keep_extra_cols = c('processed_sequence', 'processed_nt_change', 'ligation_mh_nt', 'processed_cdr3', 'cys_character', 'phe_character'))

    # Frame calculations
    adjusted_grouped = get_cdr3_frame(adjusted_grouped)

    # Subset data by columns
    cols3 = c('v_gene', 'j_gene', 'v_frame', 'j_frame', 'v_seq_len', 'j_seq_len', 'v_trim', 'j_trim', 'ligation_mh', 'frame_type', 'frame_stop', 'processed_sequence', 'processed_cdr3', 'processed_nt_change', 'ligation_mh_nt')
    final = adjusted_grouped[, ..cols3]
    return(final)
}
