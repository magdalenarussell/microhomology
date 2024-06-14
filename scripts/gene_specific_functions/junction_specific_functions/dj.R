get_oriented_whole_nucseqs <- function(whole_nucseq = get_whole_nucseqs()){
    # reorient sequence so that it is 5 -> 3 on the actual trimmed strand
    whole_nucseq[substring(gene, 4, 4) == 'J', sequence := unlist(lapply(sequence, function(x) as.character(reverseComplement(DNAString(x)))))]
    whole_nucseq[substring(gene, 4, 4) == 'J', j_gene_sequence := sequence]
    whole_nucseq[substring(gene, 4, 4) == 'D', d_gene_sequence := sequence]
    return(whole_nucseq[, -c('sequence')])
}

source(paste0(MOD_PROJECT_PATH,'/scripts/gene_count_specific_functions/double.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/mh_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/gene_specific_functions/junction_specific_functions/general_junction.R'))

get_processed_sequences <- function(adjusted_df, frame_data, pnucs = 2){
    stopifnot('ligation_mh' %in% colnames(adjusted_df))
    # d_gene_sequence should be the bottom strand oriented 5>3
    cols = colnames(adjusted_df)
    # get overlapping sequences
    overlaps = adjusted_df[, get_overlapping_regions('d_gene', 'j_gene', d_gene_sequence, j_gene_sequence, d3_trim, j_trim, ligation_mh, positions = c('up', 'mid', 'down'), rename_cols = FALSE)]
    tog = cbind(adjusted_df, overlaps) 

    # get bottom gene sequence 
    j_gene_bottom_seq = reorient_j_bottom_strand(adjusted_df$j_gene_sequence)

    # get pnucs
    require(Biostrings)
    d_pnucs = get_pnucs(adjusted_df$d_gene_sequence, 'top', pnucs)    
    j_pnucs = get_pnucs(j_gene_bottom_seq, 'bottom', pnucs)    

    d_gene_top_seq_p = paste0(adjusted_df$d_gene_sequence, d_pnucs)
    j_gene_bottom_seq_p = paste0(j_pnucs, j_gene_bottom_seq)

    j_gene_top_seq_p = as.character(complement(DNAStringSet(j_gene_bottom_seq_p)))

    tog$d_gene_sequence_pnucs = d_gene_top_seq_p
    tog$j_gene_sequence_pnucs_top = j_gene_top_seq_p

    # trim sequences
    tog[, trimmed_d_gene_sequence := substring(d_gene_sequence_pnucs, 1, nchar(d_gene_sequence_pnucs) - 2 - d3_trim)]
    tog[, trimmed_j_gene_sequence := substring(j_gene_sequence_pnucs_top, j_trim + 3)]

    # check MH
    tog[ligation_mh > 0, d_gene_mh := substring(trimmed_d_gene_sequence, nchar(trimmed_d_gene_sequence) - ligation_mh + 1)]
    tog[ligation_mh > 0, j_gene_mh := substring(trimmed_j_gene_sequence, 1, ligation_mh)]
    stopifnot(all(tog[ligation_mh > 0]$d_mid == tog[ligation_mh > 0]$d_gene_mh))
    stopifnot(all(tog[ligation_mh > 0]$d_mid == tog[ligation_mh > 0]$j_gene_mh))

    # get sequences corrected for MH
    tog[, trimmed_j_gene_sequence_wo_mh := substring(trimmed_j_gene_sequence, ligation_mh + 1)]

    # get processed sequences
    tog[, processed_sequence := paste(trimmed_d_gene_sequence, trimmed_j_gene_sequence_wo_mh, sep = '')]
    tog[, processed_nt_change := d_trim + j_trim + ligation_mh]

    cols = c(cols, 'processed_sequence', 'processed_nt_change')
    return(tog[, ..cols])
}

get_frames_data <- function(){
    # This function will only return frame data for annotations that are observable, after re-annotation to accommodate cases of maximal ligation-mh
    # Read only necessary columns and apply filter
    frames = fread('https://raw.githubusercontent.com/phbradley/conga/master/conga/tcrdist/db/combo_xcr.tsv')[organism == 'human' & chain == substring(CHAIN_TYPE, 3, 3)]

    # Combine operations to reduce redundancy
    d = frames[region == 'D', .(d_gene = id, d_frame = frame, d_seq = nucseq)]
    j = frames[region == 'J', .(j_gene = id, j_frame = frame, j_seq = nucseq)]

    # Merge operations
    d$dummy = 1
    j$dummy = 1
    gene_pairs = merge(d, j, by = 'dummy', allow.cartesian = TRUE)

    # Get all trimming sites
    trims = data.table(expand.grid(d_trim = seq(LOWER_TRIM_BOUND, 16), 
                                   j_trim = seq(LOWER_TRIM_BOUND, 18)))
    trims$dummy = 1

    # Merge genes and trims
    all = merge(gene_pairs, trims, by = 'dummy', allow.cartesian = TRUE)[, -c('dummy')]
    # Get sequence lengths and subset data to necessary columns
    all[, c('d_seq_len', 'j_seq_len') := .(nchar(d_seq), nchar(j_seq))]
    # remove impossible d trims
    all = all[(d3_trim) <= d_seq_len]

    all = get_oriented_full_sequences(all)

    dj_cols = c('d_gene', 'j_gene', 'd_frame', 'j_frame', 'd_seq_len', 'j_seq_len', 'd3_trim', 'j_trim', 'd_gene_sequence', 'j_gene_sequence')
    all = unique(all[, ..dj_cols])

    # get possible ligation mh
    dj_lig_mat = matrix(0, nrow = nrow(all), ncol = length(seq(1, 15)))

    for (overlap in seq(ncol(dj_lig_mat))){
        lig = get_possible_ligation_mh_fixed_trim(all, overlap_count = overlap, 'd_gene_sequence', 'j_gene_sequence', 'd3_trim', 'j_trim')
        dj_lig_mat[, overlap] = lig
        print(paste0('finished getting DJ ligation configurations for ', overlap, ' MH overlap'))
    }

    # get unique vals 
    unique_dj_mh_list = apply(dj_lig_mat, 1, unique)

    # combine
    adjusted_all = data.table(all, ligation_mh = unique_dj_mh_list)

    # Expand the rows
    adjusted_grouped = as.data.table(adjusted_all[, unnest(.SD, cols = c("ligation_mh"))])

    # Get oriented full sequences and group genes by common features, also subset columns again
    cols = c(dj_cols, 'ligation_mh')
    adjusted_grouped = unique(adjusted_grouped[, ..cols])
    
    # get processed sequences
    adjusted_grouped = get_processed_sequences(adjusted_grouped, frames)

    # Get stop positions and frames
    # all processed sequences have the potential to be out-of-frame!!!
    adjusted_grouped[, frame_type := 'Out']
    adjusted_grouped[, frame_stop := FALSE] 

    # Subset data by columns
    cols3 = c('d_gene', 'j_gene', 'd_frame', 'j_frame', 'd_seq_len', 'j_seq_len', 'd3_trim', 'j_trim', 'ligation_mh', 'frame_type', 'frame_stop', 'processed_sequence', 'processed_nt_change')
    final = adjusted_grouped[, ..cols3]
    return(final)
}
