get_oriented_full_sequences <- function(subject_data, whole_nucseq = get_oriented_whole_nucseqs(), gene_type = GENE_NAME){
    together = merge(subject_data, whole_nucseq, by.x = gene_type, by.y = 'gene')
    return(together)
}
