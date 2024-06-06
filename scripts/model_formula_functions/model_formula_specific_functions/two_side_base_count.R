get_all_base_variables <- function(side, trim_type = TRIM_TYPE){
    stopifnot(side %in% c('5end', '3end'))
    bases = c('GC', 'AT')
    vars = paste0(trim_type, '_', side, '_base_count_', bases)
    return(vars)
}

get_all_base_prop_length_interaction_variables <- function(side, trim_type = TRIM_TYPE){
    stopifnot(side %in% c('5end', '3end'))
    bases = c('GC', 'AT')
    length = paste0(substring(trim_type, 1, 1), '_length')
    vars = paste0(trim_type, '_', side, '_base_count_', bases, '_prop_', length, '_interaction')
    return(vars)
}

count_bases_seq_list <- function(seq_list, side, trim_type = TRIM_TYPE){
    subset = unique(seq_list)
    bases = c('A', 'T', 'G', 'C')
    registerDoParallel(cores=NCPU)  
    cl = makeCluster(NCPU, type="FORK")
    seq_list_DNA = parLapply(cl, subset, function(x) DNAString(x))
    counts = parLapply(cl, seq_list_DNA, function(x) letterFrequency(x, letters = bases))
    counts_dt = rbindlist(lapply(counts, as.data.frame.list))
    counts_dt[, 'AT' := A+T]
    counts_dt[, 'GC' := G+C]
    stopCluster(cl)  
    paired_cols = c()
    for (base in c('AT', 'GC')) {
        new_name = paste0(trim_type, '_', side, '_base_count_', base)
        setnames(counts_dt, base, new_name)
        paired_cols = c(paired_cols, new_name)
    }
    together = cbind(subset, counts_dt[, ..paired_cols])
    setnames(together, 'subset', paste0(trim_type, '_', side, '_seq'))
    return(together)
}
    
get_left_right_seq_vars <- function(motif_data, left_nuc_count = LEFT_SIDE_TERMINAL_MELT_LENGTH, beyond_motif, single_stranded = FALSE, whole_nucseq = get_oriented_whole_nucseqs(), gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    trims = seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND)
    
    gene_col = whole_nucseq$gene
    genes = gene_col[substring(gene_col, 4, 4) == toupper(substring(gene_type, 1, 1))]
    together = data.table(gene = rep(genes, each = length(trims)))
    together[[trim_type]] = rep(trims, length(genes))
    together = merge(together, whole_nucseq, by = 'gene')

    setnames(together, 'gene', gene_type)

    cols = c(trim_type, paste0(gene_type, '_sequence') , paste0(gene_type))
    together = together[, ..cols]

    if ('name' %in% colnames(together)){
        together = together[, -c('name')]
    }

    # get terminal seq
    if (isTRUE(beyond_motif)){
        right_position_shift = RIGHT_NUC_MOTIF_COUNT
        left_position_shift = LEFT_NUC_MOTIF_COUNT 
    } else {
        right_position_shift = 0
        left_position_shift = 0
    }

    require(Biostrings)
    seq = DNAStringSet(together[[paste0(gene_type, '_sequence')]])
    if (PNUC_COUNT > 0){
        possible_pnucs_5_to_3 = substring(reverseComplement(seq),1, PNUC_COUNT)
    } else if (PNUC_COUNT < 0){
        possible_pnucs_5_to_3 = DNAString() 
        seq = substring(seq, 1, nchar(seq) + PNUC_COUNT)             
    } else {
        possible_pnucs_5_to_3 = DNAString() 
    }

    together[[paste0(gene_type, '_sequence_pnuc')]] = paste(as.character(seq), as.character(possible_pnucs_5_to_3), sep = '')

    if (isFALSE(single_stranded)){
        together[, paste0(trim_type, '_3end_seq') := substring(get(paste0(gene_type, '_sequence_pnuc')), nchar(get(paste0(gene_type, '_sequence_pnuc'))) - get(trim_type) + 1 + right_position_shift - abs(PNUC_COUNT), nchar(get(paste0(gene_type, '_sequence_pnuc')))-2*abs(PNUC_COUNT))]
    } else {
        together[, paste0(trim_type, '_3end_seq') := substring(get(paste0(gene_type, '_sequence_pnuc')), nchar(get(paste0(gene_type, '_sequence_pnuc'))) - get(trim_type) + 1 + right_position_shift - abs(PNUC_COUNT), nchar(get(paste0(gene_type, '_sequence_pnuc'))))]
    }

    if (is.numeric(left_nuc_count)){
        together[, paste0(trim_type, '_5end_seq') := substring(get(paste0(gene_type, '_sequence_pnuc')), nchar(get(paste0(gene_type, '_sequence_pnuc'))) - (get(trim_type) + left_nuc_count) + 1 - abs(PNUC_COUNT), nchar(get(paste0(gene_type, '_sequence_pnuc')))-get(trim_type) - left_position_shift - abs(PNUC_COUNT))]
    } else if (left_nuc_count == 'right_nuc_count') {
        together[, paste0(trim_type, '_5end_seq') := substring(get(paste0(gene_type, '_sequence_pnuc')), nchar(get(paste0(gene_type, '_sequence_pnuc'))) - (get(trim_type) + nchar(get(paste0(trim_type, '_3end_seq')))) + 1 - abs(PNUC_COUNT), nchar(get(paste0(gene_type, '_sequence_pnuc')))-get(trim_type) - left_position_shift - abs(PNUC_COUNT))]
    }
   
    cols = colnames(together)[!(colnames(together) %like% '_sequence')]
    motif_data_together = merge(motif_data, unique(together[, ..cols]), by = c(paste0(gene_type), trim_type))
    return(motif_data_together)
}

process_for_two_side_base_count <- function(motif_data, left_nuc_count = LEFT_SIDE_TERMINAL_MELT_LENGTH, beyond_motif = FALSE, single_stranded = FALSE, whole_nucseq = get_oriented_whole_nucseqs(), gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    vars = get_all_base_variables('5end', trim_type)
    if (!all(vars %in% colnames(motif_data))){
        motif_data = get_left_right_seq_vars(motif_data, left_nuc_count, beyond_motif, single_stranded, whole_nucseq, gene_type, trim_type)

        for (side in c('5end', '3end')){
            col = paste0(trim_type, '_', side, '_seq')
            base_counts = count_bases_seq_list(motif_data[[col]], side, trim_type)
            motif_data = merge(motif_data, base_counts, by = col)
        }
    }

    cols = colnames(motif_data)[!(colnames(motif_data) %like% '_seq')]
    if ('processed_sequence' %in% colnames(motif_data)){
        cols = c(cols, 'processed_sequence')
    }
    return(motif_data[, ..cols])
}
