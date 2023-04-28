get_whole_nucseqs <- function(){
    whole_nucseq = fread(get(paste0('WHOLE_NUCSEQS_alpha')))[, -c('name')]
    return(whole_nucseq[, c('gene', 'sequence')])
}

get_subsequence <- function(seq_list, length, gene_type){
    if (toupper(gene_type) == 'J' | toupper(gene_type) == 'D5'){
        subseqs = substring(seq_list, 1, length) 
    } else if (toupper(gene_type) == 'V' | toupper(gene_type) == 'D3'){
        subseqs = sapply(seq_list, function(x) substring(x, nchar(x) - length + 1, nchar(x)))
    }
    return(subseqs)
}

get_similarity_distance <- function(nt_count, gene_type = 'J', whole_nucseqs){
    require(DECIPHER)
    require(Biostrings)
    tr = whole_nucseqs[gene %like% gene_type]
    if (nt_count == 'all') {
        tr[, short := sequence]
    } else {
        tr[, short:= get_subsequence(sequence, nt_count, gene_type)]
    }
    dists = DistanceMatrix(DNAStringSet(tr$short))

    tr[, dist := rowMeans(dists)]
    tr[, dist_nt_count := nt_count]
    return(tr)
}


get_pairwise_hamming <- function(joining_genes, nt_count, alignment, gene_type = substring(JOINING_GENE, 1, 1)){
    require(DECIPHER)
    require(Biostrings)
    seqs = get_whole_nucseqs()
    seqs = seqs[gene %in% joining_genes]
    if (alignment == 'align'){
        require(msa)
        a = msa(DNAStringSet(seqs$sequence), 'Muscle')
        seqs$sequence = as.character(a)
    }
    if (nt_count == 'all') {
        seqs[, short := sequence]
    } else {
        seqs[, short:= get_subsequence(sequence, nt_count, gene_type)]
    }
    dists = as.data.table(DistanceMatrix(DNAStringSet(seqs$short)))
    colnames(dists) = seqs$gene
    dists[[paste0(JOINING_GENE, '.x')]] = seqs$gene
    
    col = paste0(JOINING_GENE, '.x')
    #transform
    longer = dists %>%
        pivot_longer(!any_of(col), names_to = paste0(JOINING_GENE, '.y'), values_to = 'dist') %>%
        as.data.table()
    return(longer)
}

get_pairwise_msa <- function(joining_genes, nt_count, gene_type = substring(JOINING_GENE, 1, 1)){
    require(Biostrings)
    seqs = get_whole_nucseqs()
    seqs = seqs[gene %in% joining_genes]

    if (nt_count == 'all') {
        seqs[, short := sequence]
    } else {
        seqs[, short:= get_subsequence(sequence, nt_count, gene_type)]
    }
    mat = nucleotideSubstitutionMatrix(match = 1, mismatch = -2, baseOnly = TRUE)
    pairwise_aligns = data.table()
    for (ind in seq(nrow(seqs))){
        gene1 = seqs[ind]$gene
        for (ind2 in seq(nrow(seqs))){
            gene2 = seqs[ind2]$gene
            if (nrow(pairwise_aligns) > 0){
                if (nrow(pairwise_aligns[get(paste0(JOINING_GENE, '.y')) == gene1 & get(paste0(JOINING_GENE, '.x')) == gene2]) > 0){
                    next
                }
            }
            score = pairwiseAlignment(DNAStringSet(seqs$short)[ind], DNAStringSet(seqs$short)[ind2], scoreOnly = TRUE, substitutionMatrix = mat, gapOpening = 4, gapExtension = 2)
            temp = data.table(pairwise_score = score)
            temp[[paste0(JOINING_GENE, '.x')]] = gene1
            temp[[paste0(JOINING_GENE, '.y')]] = gene2
            pairwise_aligns = rbind(pairwise_aligns, temp)
        }
    }
    return(pairwise_aligns)
}

transform_cluster_data_to_pairwise <- function(cluster_data){
    cols = c(GENE_NAME, JOINING_GENE, 'cluster')
    condensed = unique(cluster_data[, ..cols])

    new_dt = data.table()
    
    for (v in unique(condensed[[GENE_NAME]])){
        temp = condensed[get(GENE_NAME) == v]
        for (clust in unique(temp$cluster)){
            clust_temp = temp[cluster == clust]
            joining_set = unique(clust_temp[[JOINING_GENE]])
            all_joining = unique(temp[[JOINING_GENE]])
            for (j1 in all_joining){
                for (j2 in all_joining){
                    if (nrow(new_dt) > 0){
                        if (nrow(new_dt[get(paste0(JOINING_GENE, '.y')) == j1 & get(paste0(JOINING_GENE, '.x')) == j2])){
                            next
                        }
                    }
                    if (j1 %in% joining_set & j2 %in% joining_set){
                        new_temp = data.table(cluster = clust)
                        new_temp[[paste0(JOINING_GENE, '.x')]] = j1
                        new_temp[[paste0(JOINING_GENE, '.y')]] = j2
                        new_temp[[GENE_NAME]] = v
                    } else {
                        new_temp = data.table(cluster = 0)
                        new_temp[[paste0(JOINING_GENE, '.x')]] = j1
                        new_temp[[paste0(JOINING_GENE, '.y')]] = j2
                        new_temp[[GENE_NAME]] = v
                    }
                    new_dt = rbind(new_dt, new_temp)
                }
            }
        }
    }
    return(new_dt)
}

