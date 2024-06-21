get_possible_pnucs <- function(data){
    require(Biostrings)
    seqs = get_whole_nucseqs()

    # get vj insert nucs
    data = get_insertion_indices(data)

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

get_v_cdr3_nucseq <- function(){
    frames = fread('https://raw.githubusercontent.com/phbradley/conga/master/conga/tcrdist/db/combo_xcr.tsv')[organism %like% 'human' & chain == CHAIN_SUBTYPE][substring(id, 3, 3) == substring(CHAIN_TYPE, 3, 3)]

    vframe = frames[id %like% 'V']
    vframe[, c("extra", "cdr1", "cdr2", 'cdr3') := tstrsplit(cdr_columns, ";", fixed = TRUE)]
    vframe[, c("cdr3_start", "CDR3_end") := tstrsplit(cdr3, "-", fixed = TRUE)]
    vframe[, seq_start := substring(aligned_protseq, 1, cdr3_start)]
    vframe[, dot_count := stringr::str_count(seq_start, "\\.")]
    vframe[, cys_start_index := as.numeric(cdr3_start) - dot_count]
    setnames(vframe, 'id', 'v_gene')
    
    vframe[, v_cdr3_nucseq := substring(nucseq, cys_start_index*3-2)]
    vframe[, v_cdr3_prot := as.character(translate(DNAStringSet(v_cdr3_nucseq)))]

    # check work
    vframe[, c("extra_prot", "cdr1_prot", "cdr2_prot", 'cdr3_prot') := tstrsplit(cdrs, ";", fixed = TRUE)]
    vframe[, cdr3_prot := stringr::str_remove_all(cdr3_prot, '\\.')]

    vframe_subset = vframe[cdr3_prot == v_cdr3_prot]

    cols = c('v_gene', 'v_cdr3_nucseq')
    return(vframe_subset[, ..cols])
}

get_insertion_indices <- function(data){
    nucseqs = get_v_cdr3_nucseq()

    data = merge(data, nucseqs)

    # get trimmed V nucseq
    data[v_trim > 0, v_cdr3_nucseq_trimmed := substring(v_cdr3_nucseq, 1, nchar(v_cdr3_nucseq) - v_trim)]
    data[v_trim == 0, v_cdr3_nucseq_trimmed := v_cdr3_nucseq]

    # get n_index
    data[, v_cdr3_nchar := nchar(v_cdr3_nucseq_trimmed)]
    data[, n_index := v_cdr3_nchar]

    # get n insertions
    data[, vj_insert_nucs := toupper(substring(cdr3_nucseq, v_cdr3_nchar, v_cdr3_nchar + vj_insert-1))]

    return(data)
}

reformat_data <- function(data){
    stopifnot(LOCUS %in% c('alpha', 'gamma'))
    if (range(data$v_trim)[1] >= 0){
        data = get_possible_pnucs(data)
    }
    setnames(data, 'd0_trim', 'd5_trim', skip_absent = TRUE)
    setnames(data, 'd1_trim', 'd3_trim', skip_absent = TRUE)
    return(data)
}
