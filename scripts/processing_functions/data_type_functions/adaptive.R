filter_chain <- function(data){
    stopifnot(LOCUS %like% 'TR')
    if (LOCUS == 'TRA'){
        rearr = 'VJ'
    } else {
        rearr = 'VDJ'
    }
    data = data[rearrangement_type == rearr]
    return(data)
}

convert_inserts <- function(data){
    stopifnot(LOCUS %like% 'TR')
    if (LOCUS == 'TRA'){
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
    if (data[[gene_col]][1] %like% 'TCR'){
        mapping = fread('https://raw.githubusercontent.com/kmayerb/tcrdist3/master/tcrdist/db/adaptive_imgt_mapping.csv')[species == 'human']
        gene_type = toupper(substring(gene_col, 1, 1))
        data[, c(gene_col, "allele") := tstrsplit(get(gene_col), "\\*0")]
        tog = merge(data, mapping, by.x = gene_col, by.y = 'adaptive', all.x = TRUE, allow.cartesian = TRUE)
        setnames(tog, 'imgt', paste0(tolower(gene_type), '_gene'))
        return(tog[, -c('allele', 'species')])
    } else {
        return(data)
    }
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
            data = convert_trims(data, paste0(trim, 'deletions'))
            data[[paste0(trim, 'trim')]] = as.numeric(data[[paste0(trim, 'trim')]])
        }

        data = convert_inserts(data)
        data = filter_chain(data)
    }
    return(data)
}

reformat_data <- function(data){
    data = convert_adaptive_style_to_imgt(data) 
    return(data)
}
