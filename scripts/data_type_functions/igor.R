adjust_pnucs <- function(data){
    if ('v_pnuc' %in% colnames(data)){
        data[v_pnuc > 0, v_trim := v_trim - v_pnuc]
    }
    if ('d0_pnuc' %in% colnames(data)){
        data[d0_pnuc > 0, d0_trim := d0_trim - d0_pnuc]
    }
    if ('d1_pnuc' %in% colnames(data)){
        data[d1_pnuc > 0, d1_trim := d1_trim - d1_pnuc]
    }
    if ('j_pnuc' %in% colnames(data)){
        data[j_pnuc > 0, j_trim := j_trim - j_pnuc]
    }
    return(data)
}

reformat_data <- function(data){
    data = adjust_pnucs(data)
    setnames(data, 'd0_trim', 'd5_trim', skip_absent = TRUE)
    setnames(data, 'd1_trim', 'd3_trim', skip_absent = TRUE)
    return(data)
}
