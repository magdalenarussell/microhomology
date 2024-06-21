reformat_data <- function(data){
    stopifnot(LOCUS == 'alpha')
    if (all(unique(data$productive) %in% c(TRUE, FALSE, NA))){
        data[, productive := as.character(productive)]
        data[productive == 'TRUE', productive := 'productive']
        data[productive == 'FALSE', productive := 'nonproductive']
    } 

    data = get_vdjserver_genes(data)
    data = get_vdjserver_trims(data)
    data = get_possible_pnucs(data)

    cols = c('v_gene', 'j_gene', 'v_trim', 'j_trim', 'vj_insert', 'productive')
    subset = data[, ..cols]
    setnames(data, 'd0_trim', 'd5_trim', skip_absent = TRUE)
    setnames(data, 'd1_trim', 'd3_trim', skip_absent = TRUE)
    return(subset)
}

get_vdjserver_genes <- function(data) {
    data[, v_gene := v_call]
    data[, j_gene := j_call]
    data[v_gene %like% ',', v_gene := str_split_i(v_gene, ',', 1)]
    data[j_gene %like% ',', j_gene := str_split_i(j_gene, ',', 1)]

    data = data[v_gene %like% 'TRAV' & j_gene %like% 'TRAJ']
    data = data[!is.na(v_gene)][!(is.na(j_gene))]
    return(data)
}

get_vdjserver_trims <- function(data){
    data[, v_trim := substring(v_cigar, nchar(v_cigar) - 3, nchar(v_cigar))]
    data[, j_trim := substring(j_cigar, 1, 6)]

    for (gt in c('v_', 'j_')){
        data[!(get(paste0(gt, 'trim')) %like% 'N'), paste0(gt, 'trim') := 0]
        data[get(paste0(gt, 'trim')) %like% 'S', paste0(gt, 'trim') := str_split_i(get(paste0(gt, 'trim')), 'S', 2)]
        data[get(paste0(gt, 'trim')) %like% 'N', paste0(gt, 'trim') := str_split_i(get(paste0(gt, 'trim')), 'N', 1)]
        data[, paste0(gt, 'trim') := as.numeric(get(paste0(gt, 'trim')))]
    }

    data[, j_pnuc := p5j_length]
    data[, v_pnuc := p3v_length]

    # get insertions
    data[, vj_insert := np1_length]
    data[, vj_insert_nucs := np1]

    data = data[!is.na(v_trim)][!is.na(j_trim)][!is.na(vj_insert)]
    return(data)
}

get_possible_pnucs <- function(data){
    require(Biostrings)
    seqs = get_whole_nucseqs()

    for (gt in c('v_', 'j_')){
        if (!all(is.na(data[[paste0(gt, 'pnuc')]]))){
            next
        }
        
        setnames(seqs, colnames(seqs)[colnames(seqs) != 'gene'], paste0(gt, 'whole_seq'))
        data = merge(data, seqs, by.x = paste0(gt, 'gene'), by.y = 'gene', all.x = TRUE)
        
        # get possible pnucs
        if (gt == 'v_'){
            data[v_trim == 0 & vj_insert > 0, v_possible_pnucs := substring(reverseComplement(DNAStringSet(v_whole_seq)),1, 2)]
        } else if (gt == 'j_'){
            data[j_trim == 0 & vj_insert > 0, j_possible_pnucs := substring(reverseComplement(DNAStringSet(j_whole_seq)), nchar(j_whole_seq) - 1, nchar(j_whole_seq))]
        }

        # see if pnucs are present
        data[!is.na(get(paste0(gt, 'possible_pnucs'))) & get(paste0(gt, 'trim')) == 0 & vj_insert > 0, paste0(gt, 'pnuc_nt') := get_pnucs_from_vj_nucs(get(paste0(gt,'possible_pnucs')), vj_insert_nucs, gt)]

        # make pnucs nonzero trims
        data[nchar(get(paste0(gt, 'pnuc_nt'))) > 0, paste0(gt, 'trim') := -1* nchar(get(paste0(gt, 'pnuc_nt')))]
    }

    data = correct_overlapping_pnucs(data) 

    data[v_trim < 0, vj_insert := vj_insert + v_trim] 
    data[j_trim < 0, vj_insert := vj_insert + j_trim] 
    return(data)
}

get_pnucs_from_vj_nucs <- function(list_of_pnucs, list_of_inserts, pnuc_type) {
    stopifnot(pnuc_type %in% c('v_', 'j_'))
        
    if (pnuc_type == 'v_'){
        pnuc_inserts = substring(list_of_inserts, 1, 2)
        first_index = 1
        second_index = 2
    } else if (pnuc_type == 'j_'){
        pnuc_inserts = substring(list_of_inserts, nchar(list_of_inserts)-1, nchar(list_of_inserts))
        first_index = 2
        second_index = 1 
    }

    pnuc = character(length(list_of_pnucs))

    # check if both characters in the same position match
    first_pnuc_match = substr(list_of_pnucs, first_index, first_index) == substr(pnuc_inserts, first_index, first_index)

    # add matching first characters to the vector
    pnuc[first_pnuc_match] = substr(list_of_pnucs[first_pnuc_match], first_index, first_index)

    # check if second characters match
    second_pnuc_match = substr(list_of_pnucs, second_index, second_index) == substr(pnuc_inserts, second_index, second_index)

    # add matching second characters to the vector
    if (first_index == 1){
        pnuc[first_pnuc_match & second_pnuc_match] = paste0(pnuc[first_pnuc_match & second_pnuc_match], substr(list_of_pnucs[first_pnuc_match & second_pnuc_match], second_index, second_index))
    } else if (first_index == 2){
        pnuc[first_pnuc_match & second_pnuc_match] = paste0(substr(list_of_pnucs[first_pnuc_match & second_pnuc_match], second_index, second_index),pnuc[first_pnuc_match & second_pnuc_match])
    }

    # add empty strings for non-matches
    pnuc[!first_pnuc_match] = ""

    return(pnuc)
}

correct_overlapping_pnucs <- function(data){
    data[v_trim == -2 & j_trim == -1 & nchar(vj_insert_nucs) == 2, v_trim := v_trim + 1]
    data[v_trim == -1 & j_trim == -2 & nchar(vj_insert_nucs) == 2, j_trim := j_trim + 1]
    data[v_trim == -1 & j_trim == -1 & nchar(vj_insert_nucs) == 1, j_trim := j_trim + 1]
    data[v_trim == -2 & j_trim == -2 & nchar(vj_insert_nucs) == 2, c("v_trim", "j_trim") := list(v_trim + 1, j_trim + 1)]
    data[v_trim == -2 & j_trim == -2 & nchar(vj_insert_nucs) == 3, j_trim := j_trim + 1]
    return(data)
}

