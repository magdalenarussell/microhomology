get_gene_sequence <- function(whole_nucseqs, gene_name, gene_seq_length, pnuc_count = 2){
    whole_nucseq = whole_nucseqs[toupper(substring(gene, 4, 4)) == toupper(substring(GENE_NAME, 1,1))]
    setnames(whole_nucseq, 'gene', GENE_NAME)
    gene = whole_nucseq[get(GENE_NAME) == gene_name]

    # get sequence
    require(Biostrings)
    if (GENE_NAME == 'v_gene'){
        whole_gene_seq = DNAString(gene$sequence)
    } else if (GENE_NAME == 'j_gene'){
        whole_gene_seq = DNAString(gene$sequence)
        whole_gene_seq = reverseComplement(whole_gene_seq)
    } 
    possible_pnucs_5_to_3 = substring(reverseComplement(whole_gene_seq),1, pnuc_count)
    whole_gene_with_pnucs = c(unlist(whole_gene_seq), unlist(possible_pnucs_5_to_3))
    subset = substring(whole_gene_with_pnucs, nchar(whole_gene_with_pnucs) - (gene_seq_length + pnuc_count-1), nchar(whole_gene_with_pnucs))
    return(subset)
}

get_plot_positions_for_gene_sequence <- function(gene_sequence, pnuc_count = 2){
    end_position = -1 * pnuc_count + 0.5 
    start_position = nchar(gene_sequence) - pnuc_count - 0.5
    positions = seq(start_position, end_position, by = -1)
    together = data.table(base = str_split(as.character(unlist(gene_sequence)), '')[[1]], position = positions)
    return(together)
}


