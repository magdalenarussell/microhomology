get_pwm_score <- function(pwm, motif, positions){
    motif_elements = strsplit(motif, "")[[1]] 
    names(motif_elements) = positions
    score = 0
    for (index in 1:length(motif_elements)){
        element = motif_elements[index]
        element_score = pwm[base == element & parameter == names(element)]$coefficient
        score = score + element_score 
    }
    return(score)
}

get_base_count_pwm_score <- function(pwm, row, trim_type = TRIM_TYPE){
    score = 0
    for (count in c(paste0(trim_type, '_left_base_count_GC'), paste0(trim_type, '_right_base_count_AT'), paste0(trim_type, '_right_base_count_GC'))){
        element = row[[count]]
        element_score = (pwm[parameter == count]$coefficient) * as.numeric(element)
        score = score + element_score
    }
    return(score)
}

get_all_base_count_pwm_score <- function(predicted_trims, pwm, trim_type = TRIM_TYPE){
    counts = c(paste0(trim_type, '_left_base_count_GC'), paste0(trim_type, '_right_base_count_AT'), paste0(trim_type, '_right_base_count_GC'))
    condensed = unique(predicted_trims[,..counts])
    condensed[, pwm_score := apply(condensed, 1, get_base_count_pwm_score, pwm = pwm, trim_type = trim_type)]
    return(condensed)
}

get_all_pwm_score <- function(predicted_trims, pwm, per_gene_resid, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    positions = get_positions(trim_type)
    cols = c(paste0(gene_type, '_group'), trim_type, paste0(trim_type, '_motif'), positions)
    condensed = unique(predicted_trims[,..cols]) 
    # calculate pwm score by gene and trim length
    condensed[, pwm_score := unlist(lapply(get(paste0(trim_type, '_motif')), function(x) as.numeric(get_pwm_score(pwm, x, positions))))]
    if (!is.null(per_gene_resid)){
        tog = merge(condensed, per_gene_resid)
    } else {
        tog = condensed
    }
    return(tog)
}

compare_weight_by_motif <- function(predicted_trims, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    cols = c(paste0(gene_type, '_group'), trim_type, paste0(trim_type, '_motif'))
    avg = predicted_trims[, sum(weighted_observation), by = cols]
    setnames(avg, 'V1', 'total_weight')
    cols2 = c(paste0(gene_type, '_group'), paste0(trim_type, '_motif'))
    avg2 = avg[, sum(total_weight), by = cols2]
    setnames(avg2, 'V1', 'total_weight_by_gene')
    return(avg2)
}
