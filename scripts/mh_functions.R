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

get_mh_ligation_mh <- function(seq1, seq2, aligning_trim){
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
    if (aligning_trim == 'j_trim' | aligning_trim == 'd5_trim'){
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

inner_mh_function <- function(X, Y){
    if (nchar(X) != nchar(Y)){
        mh = 0
    } else if (nchar(X) == 0 | nchar(Y) == 0){
        mh = 0
    } else {
        # Split strings into characters
        x_chars <- str_split(X, "", simplify = TRUE)
        y_chars <- str_split(Y, "", simplify = TRUE)
        
        # Compare characters and count matches
        mh = sum(x_chars == y_chars)
    }
    return(mh)
}

get_mh <- function(seq1, seq2){
    seq2_comp = as.character(complement(DNAStringSet(seq2)))
    mh = mcmapply(function(X,Y) inner_mh_function(X,Y), seq1, seq2_comp)
    return(mh)
}

get_fully_contiguous_mh <- function(seq1, seq2){
    noncontig_mh = get_mh(seq1, seq2)
    contig_mh = ifelse(noncontig_mh == nchar(seq1), noncontig_mh, 0)
    return(contig_mh)
}
