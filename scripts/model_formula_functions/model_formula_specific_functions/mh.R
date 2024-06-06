source(paste0(MOD_PROJECT_PATH,'/scripts/mh_functions.R'))

reorient_j_bottom_strand <- function(j_gene_5_3_bottom_strand){
    # return j_gene bottom strand (which was previously oriented 5' > 3' for consistency with the v_gene) to be 3' > 5'
    require(stringi)
    return(stri_reverse(j_gene_5_3_bottom_strand))
}

get_overlap_names <- function(overlap_count, positions){
    names = c()
    for (pos in positions){
        temp = c(paste0('v_gene_', pos, '_overlap_', overlap_count), 
              paste0('j_gene_', pos, '_overlap_', overlap_count))
        names = c(names, temp)
    }
    return(names)
}

get_all_mh_prop_variables <- function(overlap_vector, pos = c('up', 'mid', 'down')){
    vars = c()
    for (o in overlap_vector){
        for (p in pos){
            if (o == 0 & p == 'mid'){
                next
            }
            var = paste0('mh_prop_', p, '_overlap_', o) 
            vars = c(vars, var)
        }
    }
    return(vars)
}

get_all_mh_count_variables <- function(overlap_vector, pos = c('up', 'mid', 'down')){
    vars = c()
    for (o in overlap_vector){
        for (p in pos){
            if (o == 0 & p == 'mid'){
                next
            }
            var = paste0('mh_count_', p, '_overlap_', o) 
            vars = c(vars, var)
        }
    }
    return(vars)
}


get_all_mh_prop_length_interaction_variables <- function(overlap_vector){
    pos = c('up', 'down')
    lengths = c('j_length', 'v_length')
    vars = c()
    for (o in overlap_vector){
        for (index in seq(pos)){
            p = pos[index]
            len = lengths[index]
            var = paste0('mh_prop_', p, '_overlap_', o, '_', len, '_interaction') 
            vars = c(vars, var)
        }
    }
    return(vars)
}

get_all_mh_count_length_interaction_variables <- function(overlap_vector){
    pos = c('up', 'down')
    lengths = c('j_length', 'v_length')
    vars = c()
    for (o in overlap_vector){
        for (index in seq(pos)){
            p = pos[index]
            len = lengths[index]
            var = paste0('mh_count_', p, '_overlap_', o, '_', len, '_interaction') 
            vars = c(vars, var)
        }
    }
    return(vars)
}

get_mh_prop <- function(seq1, seq2){
    mh = get_mh(seq1, seq2)
    lengths = nchar(seq1) 
    mh_prop = mh/lengths
    return(mh_prop)
}

get_mh_prop_cols <- function(data, overlap_count, keep_gene_seqs = FALSE, prop = TRUE, positions = c('up', 'down', 'mid')){
    # get overlapping regions
    names = get_overlap_names(overlap_count, positions)
    if (!all(names %in% colnames(data))){
        data[, paste(names) := get_overlapping_regions(v_gene_sequence, j_gene_sequence, v_trim, j_trim, overlap_count, positions = positions)]

        # get MH
        for (pos in positions){
            n = paste0('mh_prop_', pos, '_overlap_', overlap_count)
            if (isFALSE(prop)){
                n = paste0('mh_count_', pos, '_overlap_', overlap_count)
            }

            v_seq_col = paste0('v_gene_', pos, '_overlap_', overlap_count)
            j_seq_col = paste0('j_gene_', pos, '_overlap_', overlap_count)
            cols = c(v_seq_col, j_seq_col)
            subset = unique(data[, ..cols])

            subset[, paste(n) := get_mh_prop(get(v_seq_col), get(j_seq_col))]
            if (isFALSE(prop)){
                subset[, paste(n) := get_mh(get(v_seq_col), get(j_seq_col))]
            }
            subset[is.na(get(n)), paste(n) := 0]
            data = merge(data, subset, by = cols)
        }

        if (keep_gene_seqs == FALSE) {
            remove = c('v_gene_sequence', 'j_gene_sequence')
            cols = colnames(data)[!(colnames(data) %in% remove)]
            data = data[, ..cols]
        }
    }
    return(data)
}

process_for_mh <- function(motif_data, whole_nucseq = get_oriented_whole_nucseqs(), overlap_vector = c(0, 1, 2, 3, 4), trim_type = TRIM_TYPE, gene_type = GENE_NAME, prop = TRUE, positions = c('up', 'down', 'mid')){
    motif_data = get_oriented_full_sequences(motif_data, whole_nucseq, gene_type)
    
    for (overlap in overlap_vector){
        motif_data = get_mh_prop_cols(motif_data, overlap, keep_gene_seqs = TRUE, prop = prop, positions = positions)
    }

    if ('processed_sequence' %in% colnames(motif_data)){
        cols = c('processed_sequence')
    } else {
        cols = c()
    }
    cols = c(colnames(motif_data)[!(colnames(motif_data) %like% '_seq')], cols)
    return(motif_data[, ..cols])
}

get_mh_config_count <- function(motif_data){
    stopifnot('ligation_mh' %in% colnames(motif_data))
    motif_data[ligation_mh == 0, nonzero_mh_indicator := 0]
    motif_data[ligation_mh > 0, nonzero_mh_indicator := 1]
    
    motif_data[, mh_config_count := sum(nonzero_mh_indicator), by = .(v_gene, j_gene, v_trim, j_trim)]
    return(motif_data)
}
