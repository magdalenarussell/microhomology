filter_chain <- function(data){
    stopifnot(LOCUS %in% c('alpha', 'beta', 'gamma'))
    if (LOCUS %in% c('alpha', 'gamma')){
        rearr = 'VJ'
    } else {
        rearr = 'VDJ'
    }
    if ('rearrangement_type' %in% colnames(data)){
        data = data[rearrangement_type == rearr]
    }
    return(data)
}

convert_inserts <- function(data){
    stopifnot(LOCUS %in% c('alpha', 'beta', 'gamma'))
    if (LOCUS %in% c('alpha', 'gamma')){
        data[n1_insertions == 'no data', n1_insertions := 0]
        data$vj_insert = as.numeric(data$n1_insertions)
    } else {
        data[n1_insertions == 'no data', n1_insertions := 0]
        data[n2_insertions == 'no data', n2_insertions := 0]

        data$vd_insert = as.numeric(data$n1_insertions)
        data$dj_insert = as.numeric(data$n2_insertions)
    }
    return(data)
}

convert_adaptive_gene_names_to_imgt <- function(data, gene_col){
    data = copy(data) 
    stopifnot((data[[gene_col]][1] %like% 'TCR') | (data[[gene_col]][1] == 'unknown')) 
    mapping = fread('https://raw.githubusercontent.com/kmayerb/tcrdist3/master/tcrdist/db/adaptive_imgt_mapping.csv')[species == 'human']

    # replace some by hand
    mapping[adaptive == 'TCRGV07-01', imgt := 'TRGV7*01']
    mapping[adaptive == 'TCRGV06-01', imgt := 'TRGV6*01']
    mapping[adaptive == 'TCRGVB-01', imgt := 'TRGVB*01']

    gene_type = toupper(substring(gene_col, 1, 1))
    data[, c(gene_col, "allele") := tstrsplit(get(gene_col), "\\*0")]
    tog = merge(data, mapping, by.x = gene_col, by.y = 'adaptive', all.x = TRUE, allow.cartesian = TRUE)
    setnames(tog, 'imgt', paste0(tolower(gene_type), '_gene'))
    return(tog[, -c('allele', 'species')])
}

convert_frame_type_to_productivity <- function(data){
    data[frame_type %in% c('Out', 'Stop'), productive := 'nonproductive']
    data[frame_type == 'In', productive := 'productive']
    return(data[, -c('frame_type')])
}

convert_trims <- function(data, trim_col){
    trim_type = substring(trim_col, 1, 1)
    if (trim_type == 'd'){
        trim_type = substring(trim_col, 1, 2)
    }
    data[get(trim_col) == 'no data', paste(trim_col) := 0]
    setnames(data, trim_col, paste0(trim_type, '_trim'))
    return(data)
}

convert_adaptive_style_to_imgt <- function(data){
    adaptive_cols = c('v_resolved', 'd_resolved', 'j_resolved')
    if (any(adaptive_cols %in% colnames(data))){
        # convert genes 
        for (col in adaptive_cols){
            data = convert_adaptive_gene_names_to_imgt(data, col)
        } 

        # convert productivity
        data = convert_frame_type_to_productivity(data)

        # convert trims
        for (trim in c('v_', 'd3_', 'd5_', 'j_')){
            if (paste0(trim, 'deletions') %in% colnames(data)){
                data = convert_trims(data, paste0(trim, 'deletions'))
                data[[paste0(trim, 'trim')]] = as.numeric(data[[paste0(trim, 'trim')]])
            }
        }

        data = convert_inserts(data)
        data = filter_chain(data)
    }
    # filter for rearrangements that have an inferred V and J gene
    data = data[!(is.na(v_gene)) & !(is.na(j_gene))] 
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

        # get vj insert nucs
        data = get_vj_insert_nucs(data)

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

get_vj_insert_nucs <- function(data){
    data[, corrected_vj_index := n2_index + 1]
    data[corrected_vj_index == 0, corrected_vj_index := NA]
    data[, vj_index_end := corrected_vj_index + vj_insert - 1]
    data[, vj_insert_nucs := substring(rearrangement, corrected_vj_index, vj_index_end)]
    return(data)
}

correct_overlapping_pnucs <- function(data){
    data[v_trim == -2 & j_trim == -1 & nchar(vj_insert_nucs) == 2, v_trim := v_trim + 1]
    data[v_trim == -1 & j_trim == -2 & nchar(vj_insert_nucs) == 2, j_trim := j_trim + 1]
    data[v_trim == -1 & j_trim == -1 & nchar(vj_insert_nucs) == 1, j_trim := j_trim + 1]
    data[v_trim == -2 & j_trim == -2 & nchar(vj_insert_nucs) == 2, c("v_trim", "j_trim") := list(v_trim + 1, j_trim + 1)]
    data[v_trim == -2 & j_trim == -2 & nchar(vj_insert_nucs) == 3, j_trim := j_trim + 1]
    return(data)
}

convert_frequency_to_count <- function(data, convert = FALSE){
    if (convert == TRUE){
        stopifnot('frequency' %in% colnames(data))
        data[, templates := 100 * frequency]
    } else {
        data[templates == 'na', templates := NA]
        data[!is.na(templates), templates := as.numeric(templates)]
        data[templates == -1, templates := 0]
    }
    return(data)
}

reformat_data <- function(data){
    stopifnot(LOCUS %in% c('alpha', 'gamma'))
    data = convert_adaptive_style_to_imgt(data) 
    data = get_possible_pnucs(data)
    data = convert_frequency_to_count(data, convert = FALSE)
    setnames(data, 'sample_name', 'subject')
    setnames(data, 'templates', 'count')
    setnames(data, 'd0_trim', 'd5_trim', skip_absent = TRUE)
    setnames(data, 'd1_trim', 'd3_trim', skip_absent = TRUE)
    return(data)
}
