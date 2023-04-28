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

get_smoothed_pval <- function(data, xvar, yvar, facet_var){
    form = as.formula(paste0(yvar, '~', xvar))
    results = data.table()
    for (v in unique(data[[facet_var]])){
        subset = data[get(facet_var) == v]
        if (nrow(subset) <= 2){
            next
        }
        reg = lm(form, data = subset)
        result = data.table(gene = v, slope = signif(summary(reg)$coefficients[xvar, 'Estimate'], 4), pvalue = signif(summary(reg)$coefficients[xvar, 'Pr(>|t|)'], 4)) 
        results = rbind(results, result)
    }
    setnames(results, 'gene', facet_var, skip_absent = TRUE)
    return(results)
}

plot_general_scatter <- function(data, xvar, yvar, xtitle, ytitle, title, facet_var, facet_col, add_trend = FALSE){
    if (isTRUE(add_trend)){
        sig_data = get_smoothed_pval(data, xvar, yvar, facet_var)
        levs = unique(sig_data[order(pvalue)][[facet_var]])
        sig_data[[facet_var]] = factor(sig_data[[facet_var]], levels = levs)
        data[[facet_var]] = factor(data[[facet_var]], levels = levs)
        data = data[!is.na(get(facet_var))]
    }

    temp_plot = ggplot() +
        geom_point(data = data, aes(x = get(xvar), y = get(yvar)), size = 4, alpha = 0.3) +
        ggtitle(title) +
        xlab(xtitle) +
        ylab(ytitle) +
        theme_cowplot(font_family = 'Arial') + 
    theme(text = element_text(size = 40), axis.text.x=element_text(size = 30), axis.text.y = element_text(size = 30), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5), legend.position = 'bottom', legend.direction = 'horizontal', legend.justification="center") + 
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5)

    if (!is.null(facet_var)){
        temp_plot = temp_plot +facet_wrap(~get(facet_var), ncol = facet_col)
    }
    
    if (isTRUE(add_trend)){
        temp_plot = temp_plot +
            geom_smooth(method = 'lm', data = data, aes(x = get(xvar), y = get(yvar)), size = 3, se = FALSE)+
            geom_text(data = sig_data, aes(label = paste0('effect size = ', slope, '\np-value = ', pvalue)), x = Inf, y = Inf, color = 'blue', size = 8, vjust = 1.3, hjust = 1.1, nudge_x = 0.05, nudge_y = 0.05)
    }

    return(temp_plot)
}

plot_general_boxplot <- function(data, xvar, yvar, xtitle, ytitle, title, facet_var, facet_col, test_to_all = TRUE){
    require(ggpubr)
    temp_plot = ggplot(data = data, aes(x = as.factor(get(xvar)), y = get(yvar))) +
        geom_boxplot(size = 2) +
        geom_jitter(size = 3, alpha = 0.5, width = 0.1) +
        facet_wrap(~get(facet_var), ncol = facet_col) +
        ggtitle(title) +
        xlab(xtitle) +
        ylab(ytitle) +
        theme_cowplot(font_family = 'Arial') + 
    theme(text = element_text(size = 40), axis.text.x=element_text(size = 30), axis.text.y = element_text(size = 30), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5), legend.position = 'none') + 
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5)

    if (isTRUE(test_to_all)){
        means = data[, mean(get(yvar)), by = facet_var]
        temp_plot = temp_plot +
            geom_hline(data = means, aes(yintercept = V1), linetype = 2, size = 2, color = 'blue')+
            stat_compare_means(label = "p", method = "t.test", ref.group = ".all.", size = 8)
    }
    
    return(temp_plot)
}

get_xaxis_multinom <- function(){
    if (TRIM_TYPE %like% 'trim'){
        label = 'Number of trimmed nucleotides'
    } else if (TRIM_TYPE == 'zero_insert'){
        label = 'Frequency of zero insertions'
    } else if (TRIM_TYPE == 'zero_process'){
        label = 'Frequency of zero trimmed or inserted nucleotides'
    }
    return(label)
}

