get_terminal_sequences <- function(combine_by_terminal, gene_type = GENE_NAME){
    # get genes and sequences
    whole_nucseq = get_oriented_whole_nucseqs()
    seqs = whole_nucseq[substring(gene, 4, 4) == toupper(substring(gene_type, 1, 1))]
    
    map = get_common_genes_from_seqs(seqs, gene_type = gene_type)
    together = merge(seqs, map, by = gene_type)

    # get terminal sequences
    terminal_length = UPPER_TRIM_BOUND + REQUIRED_COMMON_NUCS_5 
    together[, paste0(gene_type, '_terminal_seq') := substring(get(paste0(gene_type, '_sequence')), nchar(get(paste0(gene_type, '_sequence')))- (terminal_length - 1), nchar(get(paste0(gene_type, '_sequence'))))]
    if (isTRUE(combine_by_terminal)){
        cols = c(paste0(gene_type, '_terminal_seq'), paste0(gene_type, '_group'))
        together = unique(together[, ..cols])
    }
    return(together)
}

get_distances <- function(sequences, gene_var, combine_by_terminal = TRUE, full_sequence = FALSE, align = FALSE, gene_type = GENE_NAME){
    require(ape)
    require(DECIPHER)
    if (isTRUE(full_sequence)){
        seq_list = DNAStringSet(sequences[[paste0(gene_type, '_sequence')]])
        stopifnot(isFALSE(combine_by_terminal))
        stopifnot(isTRUE(align))
        seq_list = AlignSeqs(seq_list)
    } else {
        seq_list = DNAStringSet(sequences[[paste0(gene_type, '_terminal_seq')]])
        stopifnot(isFALSE(align))
    }

    dists = DistanceMatrix(seq_list)

    colnames(dists) = sequences[[gene_var]]
    row.names(dists) = sequences[[gene_var]]
    return(dists)
}

get_gene_families <- function(cluster_count, combine_by_terminal = TRUE, full_sequence = FALSE, align = FALSE, gene_type = GENE_NAME){
    require(ape)
    require(DECIPHER)

    if (isTRUE(combine_by_terminal)){
        gene_var =paste0(gene_type, '_group')  
    } else {
        gene_var = gene_type 
    }
 
    seqs = get_terminal_sequences(combine_by_terminal, gene_type = gene_type)
    dists = get_distances(seqs, gene_var, combine_by_terminal, full_sequence, align, gene_type = gene_type)
    dist_format = as.dist(dists) 
    
    clusters = hclust(dist_format)
    clusters_grouped = cutree(clusters, k = cluster_count)
    clusters_grouped_df = as.data.frame(clusters_grouped)
    clusters_grouped_df[[gene_var]] = row.names(clusters_grouped_df)
    
    together = merge(seqs, as.data.table(clusters_grouped_df), by = gene_var)
    return(list(cluster_data = together, tree = as.phylo(clusters), dists = dist_format))
}

generate_hold_out_sample <- function(motif_data, cluster_genes, gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    var = paste0(gene_type, '_group')
    motif_data_subset = motif_data[!(get(var) %in% cluster_genes)]
    sample_data = motif_data[get(var) %in% cluster_genes]
    # updating the total_tcr, p_gene, gene_weight_type, and weighted_observation variables for the newly sampled datasets
    source(paste0(MOD_PROJECT_PATH, '/scripts/sampling_procedure_functions/', GENE_WEIGHT_TYPE, '.R'), local = TRUE)
    motif_data_subset = calculate_subject_gene_weight(motif_data_subset)
    source(paste0(MOD_PROJECT_PATH, '/scripts/sampling_procedure_functions/p_gene_pooled.R'), local = TRUE)
    sample_data = calculate_subject_gene_weight(sample_data)
    return(list(sample = sample_data, motif_data_subset = motif_data_subset))
}


