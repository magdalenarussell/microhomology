# calculate root mean square error
calculate_rmse_by_gene <- function(predicted_trims){
    total_trims = unique(predicted_trims$trim_length)
    # calculate subjects per gene
    subj_counts = predicted_trims[, .N/length(total_trims), by = .(gene)]
    setnames(subj_counts, 'V1', 'subject_count')
    
    # calculate residuals for each subject, gene, trim
    predicted_trims[, residual := empirical_prob - predicted_prob]
    predicted_trims[, residual_sq := (residual)^2]

    # calculate residual sum by subject, gene
    resid_sq_sums = predicted_trims[, sum(residual_sq), by = .(gene, subject, p_gene)]
    setnames(resid_sq_sums, 'V1', 'resid_sq_sum')
    together = merge(resid_sq_sums, subj_counts, by = 'gene')

    # calculate mean residual sum by gene, calculate root mean square error by dividing by subject count
    mean_resid_sq_sums = together[, sum(resid_sq_sum), by = .(gene, subject_count, p_gene)]
    setnames(mean_resid_sq_sums, 'V1', 'mean_resid_sq_sum')
    mean_resid_sq_sums[, rmse := sqrt(mean_resid_sq_sum/subject_count)]

    return(mean_resid_sq_sums[, -c('mean_resid_sq_sum')])
}

calculate_rmse_by_gene_trim <- function(predicted_trims){
    total_trims = unique(predicted_trims$trim_length)
    # calculate subjects per gene
    subj_counts = predicted_trims[, .N/length(total_trims), by = .(gene)]
    setnames(subj_counts, 'V1', 'subject_count')
    
    # calculate residuals for each subject, gene, trim
    predicted_trims[, residual := empirical_prob - predicted_prob]
    predicted_trims[, residual_sq := (residual)^2]
    together = merge(predicted_trims, subj_counts, by = 'gene')

    # calculate mean residual sum by gene, calculate root mean square error by dividing by subject count
    mean_resid_sq_sums = together[, sum(residual_sq), by = .(gene, trim_length, subject_count, p_gene)]
    setnames(mean_resid_sq_sums, 'V1', 'mean_resid_sq_sum')
    mean_resid_sq_sums[, rmse_by_trim := sqrt(mean_resid_sq_sum/subject_count)]

    return(mean_resid_sq_sums[, -c('mean_resid_sq_sum')])
}
