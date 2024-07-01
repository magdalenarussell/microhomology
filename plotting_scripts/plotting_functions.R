source(paste0(MOD_PROJECT_PATH, '/config/file_paths.R'))

set_color_palette <- function(model_type_list, with_params = FALSE){
    require(RColorBrewer)
    if (isTRUE(with_params)){
        model_types = model_type_list[!(model_type_list %in% c('null (0 params)', '2x4motif (18 params)'))]
    } else {
        model_types = model_type_list[!(model_type_list %in% c('null', '2x4motif'))]
    }
    colors = c(brewer.pal(7, 'Dark2'), brewer.pal(8, 'Set1'), brewer.pal(7, 'Set2'))
    colors = colors[!(colors %in% c("#E41A1C", "#FFFF33", "#FFD92F"))]
    names(colors) = model_types
    temp = c("#E41A1C", "#666666", "")
    if (isTRUE(with_params)){
        names(temp) = c('2x4motif (18 params)', 'null (0 params)')
    } else {
        names(temp) = c('2x4motif', 'null')
    }
    colors = c(colors, temp)
    return(colors)
}

get_gene_sequence <- function(gene_name, gene_seq_length, pnuc_count = 2, gene_type = GENE_NAME){
    whole_nucseq = get_oriented_whole_nucseqs()
    temp_data = whole_nucseq[toupper(substring(gene, 4, 4)) == toupper(substring(gene_type, 1,1))]
    setnames(whole_nucseq, 'gene', gene_type)
    colnames(temp_data) = c(gene_type, paste0(gene_type, '_sequence'))
    together = unique(temp_data[, c('gene', paste0(gene_type, '_sequence'))])
    gene = together[gene == gene_type][1]

    # get sequence
    whole_gene_seq = DNAString(gene[[paste0(gene_type, '_sequence')]])
    possible_pnucs_5_to_3 = substring(reverseComplement(whole_gene_seq),1, pnuc_count)
    whole_gene_with_pnucs = c(unlist(whole_gene_seq), unlist(possible_pnucs_5_to_3))
    subset = substring(whole_gene_with_pnucs, nchar(whole_gene_with_pnucs) - (gene_seq_length + pnuc_count-1), nchar(whole_gene_with_pnucs))
    return(subset)
}

get_plot_positions_for_gene_sequence <- function(gene_sequence, pnuc_count = 2){
    end_position = -1 * pnuc_count + 0.5 
    start_position = nchar(gene_sequence) - pnuc_count - 0.5
    positions = seq(start_position, end_position, by = -1)
    together = data.table(base = str_split(unlist(as.character(gene_sequence)), '')[[1]], position = positions)
    return(together)
}

get_motif_colors <- function(gene_seq_positions, motif, highlight_color){
    stopifnot(nchar(motif) == 3)
    gene_seq_positions$index = seq(1, nrow(gene_seq_positions))
    gene_seq_positions[index < 15, start_motif := base]
    gene_seq_positions[index < 15, start_motif := paste0(start_motif, gene_seq_positions$base[index + 1], gene_seq_positions$base[index + 2])] 
    for (mot in motif){
        gene_seq_positions[start_motif == mot, color := highlight_color]
        start_index = gene_seq_positions[start_motif == mot]$index
        gene_seq_positions[index == start_index + 1 | index == start_index + 2, color := highlight_color]
    }
    gene_seq_positions[color != highlight_color | is.na(color), color := 'black']
    return(gene_seq_positions)
}

plot_predicted_paired_trimming_dists_single_group <- function(data, row_var, ylim = NULL, color = 'blue', seq_text = 9){
    stopifnot(row_var %in% c('v_trim', 'j_trim'))
    other_var = c('v_trim', 'j_trim')[c('v_trim', 'j_trim') != row_var]
    important_cols = c('v_trim', 'j_trim', 'predicted_prob', 'empirical_prob', 'gene_pair', 'p_gene_pair')

    data[, gene_count := sum(count), by = gene_pair]
    data[, empirical_prob := count/gene_count]

    predicted_data = data[, ..important_cols]

    predicted_data$gene_pair = factor(predicted_data$gene_pair, levels = unique(predicted_data[order(-p_gene_pair)]$gene_pair))

    predicted_data[, long_row_var := paste0(toupper(substring(row_var, 1, 1)), '-trim = ', get(row_var))]
    predicted_data$long_row_var = factor(predicted_data$long_row_var, levels = unique(predicted_data[order(get(row_var))]$long_row_var))

    # get gene sequence
    # gene_seq = get_gene_sequence(gene_name, max(data$trim_length))
    # gene_seq_with_positions = get_plot_positions_for_gene_sequence(gene_seq)
    # gene_seq_with_positions =get_motif_colors(gene_seq_with_positions, motif_highlight, 'motif')
        
    if (!is.null(ylim)){
        max_prob = ylim
    } else {
        max_prob = max(max(predicted_data$empirical_prob), max(predicted_data$predicted_prob))
    }

    plot = ggplot(predicted_data) +
        facet_grid(cols = vars(gene_pair), rows = vars(long_row_var)) +
        geom_line(aes(x = get(other_var), y = empirical_prob), size = 1.75, alpha = 0.8, color = "gray60") +
        geom_line(aes(x = get(other_var), y = predicted_prob), size = 1.75, color = color) +
        geom_vline(xintercept = 0, color = 'black', size = 2) +
        # geom_text(data = gene_seq_with_positions, y = max_prob, aes(x = position, label = base, color = color), size = seq_text) +
        geom_text(y = max_prob-0.01, aes(x = -2.1), label = '3\'- ', size = seq_text - 1) +
        geom_text(y = max_prob-0.01, aes(x = UPPER_TRIM_BOUND + 0.1), label = ' -5\'', size = seq_text - 1) +
        xlab(paste0('Number of trimmed ', toupper(substring(other_var, 1, 1)), '-gene nucleotides')) +
        ylab('Probability') +
        theme_cowplot(font_family = 'Arial') + 
        theme(legend.position = "none", text = element_text(size = 35), axis.text.x=element_text(size = 25), axis.text.y = element_text(size = 25), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5), strip.text.y = element_text(angle = 0)) + 
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) 
        # scale_color_manual(values = c(motif = motif_highlight_color, black = 'black'))
    return(plot)
}

map_positions_to_values <- function(pwm){
    pwm[side == '5end', position_value := -1 * as.numeric(substring(position, nchar(position), nchar(position)))]
    pwm[side == '3end', position_value := as.numeric(substring(position, nchar(position), nchar(position)))]
    return(pwm)
}

plot_base_count_coefficient_heatmap_single_group <- function(model_coef_matrix, with_values = FALSE, limits = NULL, prop = FALSE){
    model_coef_matrix = model_coef_matrix[(coefficient %like% 'base_count')]

    # convert to log_10
    model_coef_matrix$log_10_pdel = model_coef_matrix$value/log(10)
    
    # order variables
    unique_bases = c('AT', 'GC')
    left_bases = unique(model_coef_matrix[side %like% '5end']$base)
    right_bases = unique(model_coef_matrix[side %like% '3end']$base)
    trim_type = unique(model_coef_matrix$trim_type)

    if (length(left_bases) < 2){
        missing = unique_bases[!(unique_bases == left_bases)]
        fake_row = data.table(value = rep(NA, length(missing)), coefficient = rep("base_count", length(missing)), base = missing, position = rep('', length(missing)), side = rep('5end', length(missing)), log_10_pdel = rep(NA,length(missing)), trim_type = trim_type)
        model_coef_matrix = rbind(model_coef_matrix, fake_row, fill = TRUE)
    } 
    if (length(right_bases) < 2){
        missing = unique_bases[!(unique_bases == right_bases)]
        fake_row = data.table(value = rep(NA, length(missing)), coefficient = rep("base_count", length(missing)), base = missing, position = rep('', length(missing)), side = rep('3end', length(missing)), log_10_pdel = rep(NA,length(missing)), trim_type = trim_type)
        model_coef_matrix = rbind(model_coef_matrix, fake_row, fill = TRUE)
    } 
    
    extended_data = model_coef_matrix

    extended_data[, long_name := paste0(side, ' base count')]
    vars = c('5end base count', '3end base count')
    if (isTRUE(prop)){
        vars = c('5end base proportion', '3end base proportion')
        extended_data[, long_name := paste0(side, ' base proportion')]
    }

    extended_data$long_name = factor(extended_data$long_name, levels = vars)

    if (is.null(limits)){
        max_val = max(abs(model_coef_matrix$log_10_pdel))
        limits = c(-max_val, max_val)
    }
    
    plot = ggplot(extended_data, aes(x=long_name, y=base, fill=log_10_pdel)) +
        geom_tile() +
        theme_cowplot(font_family = 'Arial') + 
        xlab('') +
        ylab('Base type') +
        geom_vline(xintercept = 1 + 0.5, size = 3.5, color = 'black') +
        scale_fill_distiller(palette = 'PuOr', name = 'log10(probability of deletion)', limits = limits, na.value = 'gray80') + 
        theme(text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 20), legend.position = 'bottom', legend.direction = 'horizontal', legend.justification="center") +
        guides(fill = guide_colourbar(barwidth = 35, barheight = 2))
       
    if (with_values == TRUE){
        plot = plot +
            geom_text(data = model_coef_matrix, aes(x = long_name, y = base, label = round(log_10_pdel, 3)))
    }

    return(plot)
}

plot_length_coefficient_heatmap_single_group <- function(model_coef_matrix, with_values = FALSE, limits = NULL){
    model_coef_matrix = model_coef_matrix[(coefficient %like% 'length')]
    model_coef_matrix = model_coef_matrix[!(coefficient %like% 'interaction')]

    # convert to log_10
    model_coef_matrix$log_10_pdel = model_coef_matrix$value/log(10)

    if (is.null(limits)){
        max_val = max(abs(model_coef_matrix$log_10_pdel))
        limits = c(-max_val, max_val)
    }
    
    model_coef_matrix$coefficient = 'Trimming length'

    plot = ggplot(model_coef_matrix, aes(x=coefficient, y=1, fill=log_10_pdel)) +
        geom_tile() +
        theme_cowplot(font_family = 'Arial') + 
        xlab('') +
        ylab('') +
        geom_vline(xintercept = 0.5, size = 3.5, color = 'black') +
        scale_fill_distiller(palette = 'PuOr', name = 'log10(probability of deletion)', limits = limits, na.value = 'gray80') + 
        theme(text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 20), legend.position = 'bottom', legend.direction = 'horizontal', legend.justification="center", strip.text.x = element_text(angle = 90), axis.text.y = element_blank()) +
        guides(fill = guide_colourbar(barwidth = 35, barheight = 2))
       
    if (with_values == TRUE){
        plot = plot +
            geom_text(data = model_coef_matrix, aes(x = long_name, y = base, label = round(log_10_pdel, 3)))
    }

    return(plot)
}

plot_igor_coefficient_heatmap_single_group <- function(model_coef_matrix, with_values = FALSE, limits = NULL){
    model_coef_matrix = model_coef_matrix[(coefficient %like% 'trim_prob')]

    # convert to log_10
    model_coef_matrix$log_10_pdel = model_coef_matrix$value/log(10)

    if (is.null(limits)){
        max_val = max(abs(model_coef_matrix$log_10_pdel))
        limits = c(-max_val, max_val)
    }
    
    model_coef_matrix$coefficient = 'IGoR trimming parameter'

    plot = ggplot(model_coef_matrix, aes(x=coefficient, y=1, fill=log_10_pdel)) +
        geom_tile() +
        theme_cowplot(font_family = 'Arial') + 
        xlab('') +
        ylab('') +
        scale_fill_distiller(palette = 'PuOr', name = 'log10(probability of deletion)', limits = limits, na.value = 'gray80') + 
        theme(text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 20), legend.position = 'bottom', legend.direction = 'horizontal', legend.justification="center", strip.text.x = element_text(angle = 90), axis.text.y = element_blank()) +
        guides(fill = guide_colourbar(barwidth = 35, barheight = 2))
       
    if (with_values == TRUE){
        plot = plot +
            geom_text(data = model_coef_matrix, aes(x = coefficient, y = 1, label = round(log_10_pdel, 3)), size = 10)
    }

    return(plot)
}

plot_average_mh_coefficient_heatmap_single_group <- function(model_coef_matrix, with_values = FALSE, limits = NULL){
    model_coef_matrix = model_coef_matrix[(coefficient %like% 'average_mh') | (coefficient %like% 'average_interior_mh')]

    # convert to log_10
    model_coef_matrix$log_10_pdel = model_coef_matrix$value/log(10)

    if (is.null(limits)){
        max_val = max(abs(model_coef_matrix$log_10_pdel))
        limits = c(-max_val, max_val)
    }
    
    model_coef_matrix$coefficient = 'Average number of MH nucleotides'

    plot = ggplot(model_coef_matrix, aes(x=coefficient, y=1, fill=log_10_pdel)) +
        geom_tile() +
        theme_cowplot(font_family = 'Arial') + 
        xlab('') +
        ylab('') +
        scale_fill_distiller(palette = 'PuOr', name = 'log10(probability of deletion)', limits = limits, na.value = 'gray80') + 
        theme(text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 20), legend.position = 'bottom', legend.direction = 'horizontal', legend.justification="center", strip.text.x = element_text(angle = 90), axis.text.y = element_blank()) +
        guides(fill = guide_colourbar(barwidth = 35, barheight = 2))
       
    if (with_values == TRUE){
        plot = plot +
            geom_text(data = model_coef_matrix, aes(x = coefficient, y = 1, label = round(log_10_pdel, 3)), size = 10)
    }

    return(plot)
}


plot_mh_config_count_coefficient_heatmap_single_group <- function(model_coef_matrix, with_values = FALSE, limits = NULL){
    model_coef_matrix = model_coef_matrix[(coefficient %like% 'mh_config_count')]

    # convert to log_10
    model_coef_matrix$log_10_pdel = model_coef_matrix$value/log(10)

    if (is.null(limits)){
        max_val = max(abs(model_coef_matrix$log_10_pdel))
        limits = c(-max_val, max_val)
    }
    
    model_coef_matrix$coefficient = 'MH ligation configuration count'

    plot = ggplot(model_coef_matrix, aes(x=coefficient, y=1, fill=log_10_pdel)) +
        geom_tile() +
        theme_cowplot(font_family = 'Arial') + 
        xlab('') +
        ylab('') +
        scale_fill_distiller(palette = 'PuOr', name = 'log10(probability of deletion)', limits = limits, na.value = 'gray80') + 
        theme(text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 20), legend.position = 'bottom', legend.direction = 'horizontal', legend.justification="center", strip.text.x = element_text(angle = 90), axis.text.y = element_blank()) +
        guides(fill = guide_colourbar(barwidth = 35, barheight = 2))
       
    if (with_values == TRUE){
        plot = plot +
            geom_text(data = model_coef_matrix, aes(x = coefficient, y = 1, label = round(log_10_pdel, 3)), size = 10)
    }

    return(plot)
}



plot_ligation_mh_coefficient_heatmap_single_group <- function(model_coef_matrix, with_values = FALSE, limits = NULL){
    model_coef_matrix = model_coef_matrix[(coefficient %like% 'ligation_mh')]

    # convert to log_10
    model_coef_matrix$log_10_pdel = model_coef_matrix$value/log(10)

    if (is.null(limits)){
        max_val = max(abs(model_coef_matrix$log_10_pdel))
        limits = c(-max_val, max_val)
    }
    
    model_coef_matrix$coefficient = 'Ligation MH count'

    plot = ggplot(model_coef_matrix, aes(x=coefficient, y=1, fill=log_10_pdel)) +
        geom_tile() +
        theme_cowplot(font_family = 'Arial') + 
        xlab('') +
        ylab('') +
        scale_fill_distiller(palette = 'PuOr', name = 'log10(probability of deletion)', limits = limits, na.value = 'gray80') + 
        theme(text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 20), legend.position = 'bottom', legend.direction = 'horizontal', legend.justification="center", strip.text.x = element_text(angle = 90), axis.text.y = element_blank()) +
        guides(fill = guide_colourbar(barwidth = 35, barheight = 2))
       
    if (with_values == TRUE){
        plot = plot +
            geom_text(data = model_coef_matrix, aes(x = coefficient, y = 1, label = round(log_10_pdel, 3)), size = 10)
    }

    return(plot)
}


plot_mh_coefficient_heatmap_single_group <- function(model_coef_matrix, with_values = FALSE, limits = NULL, interaction_coefs = FALSE, prop = TRUE, positions = c('up', 'mid', 'down')){
    if (isTRUE(prop)){
        mh_var = 'mh_prop'
    } else {
        mh_var = 'mh_count'
    }

    model_coef_matrix = model_coef_matrix[(coefficient %like% mh_var)]

    if (isTRUE(interaction_coefs)){
        model_coef_matrix = model_coef_matrix[(coefficient %like% 'interaction')]
    } else {
        model_coef_matrix = model_coef_matrix[!(coefficient %like% 'interaction')]
    }
    # convert to log_10
    model_coef_matrix$log_10_pdel = model_coef_matrix$value/log(10)
    
    # fill in missing
    if (all(c('up', 'down') %in% positions)){
        fake_row = data.table(value = NA, coefficient = mh_var, base = '', position = 'overlap0', side = 'mid', trim_type = '', log_10_pdel = NA)
        if (isTRUE(interaction_coefs)){
            fake_row = data.table(value = NA, coefficient = paste0(mh_var, "_length_interaction"), base = '', position = 'overlap0', side = 'mid', trim_type = '', log_10_pdel = NA)
        }
        extended_data = rbind(model_coef_matrix, fake_row, fill = TRUE)
    } else {
        extended_data = model_coef_matrix
    }

    # order variables
    levs = c()
    if ('up' %in% positions){
        extended_data[side == 'up', side_long := paste0(side, 'stream')]
        levs = c(levs, 'upstream')
    }
    if ('mid' %in% positions){
        extended_data[side == 'mid', side_long := paste0(side, 'dle')]
        levs = c(levs, 'middle')
    }
    if ('down' %in% positions){
        extended_data[side == 'down', side_long := paste0(side, 'stream')]
        levs = c(levs, 'downstream')
    }   

    extended_data$side_long = factor(extended_data$side_long, levels = levs)
    extended_data[, position_value := as.numeric(substring(position, nchar(position), nchar(position)))]

    if (is.null(limits)){
        max_val = max(abs(model_coef_matrix$log_10_pdel))
        limits = c(-max_val, max_val)
    }
    
    plot = ggplot(extended_data, aes(x=side_long, y=coefficient, fill=log_10_pdel)) +
        facet_grid(rows = vars(position_value), switch = 'y') +
        geom_tile() +
        theme_cowplot(font_family = 'Arial') + 
        xlab('Relative position of MH proportion') +
        ylab('Overlapping nucleotide count\n') +
        geom_vline(xintercept = 1 + 0.5, size = 3.5, color = 'black') +
        geom_vline(xintercept = 2 + 0.5, size = 3.5, color = 'black') +
        scale_fill_distiller(palette = 'PuOr', name = 'log10(probability of deletion)', limits = limits, na.value = 'gray80') + 
        theme(text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 20), legend.position = 'bottom', legend.direction = 'horizontal', legend.justification="center", strip.text.x = element_text(angle = 90), axis.text.y = element_blank()) +
        guides(fill = guide_colourbar(barwidth = 35, barheight = 2))
       
    if (with_values == TRUE){
        plot = plot +
            geom_text(data = model_coef_matrix, aes(x = long_name, y = base, label = round(log_10_pdel, 3)))
    }

    if (isTRUE(interaction_coefs)){
        plot = plot + 
            xlab('Relative position of MH proportion\ntrimming length interaction')
    }
    if (isFALSE(prop)){
        plot = plot + 
            xlab('Relative position of MH count')
    }

    if (length(positions) == 1){
        plot = plot + geom_vline(xintercept = 0.5, size = 3.5, color = 'black') 
    }

    return(plot)
}



plot_motif_coefficient_heatmap_single_group <- function(model_coef_matrix, with_values = FALSE, limits = NULL){
    model_coef_matrix = model_coef_matrix[coefficient == 'motif']
    model_coef_matrix = map_positions_to_values(model_coef_matrix)

    # convert to log_10
    model_coef_matrix$log_10_pdel = model_coef_matrix$value/log(10)
    # order variables
    model_coef_matrix$base = factor(model_coef_matrix$base, levels = c('T', 'G', 'C', 'A'))
    model_coef_matrix$position_value = factor(model_coef_matrix$position_value)

    if (is.null(limits)){
        max_val = max(abs(model_coef_matrix$log_10_pdel))
        limits = c(-max_val, max_val)
    }
    
    motif_length = RIGHT_NUC_MOTIF_COUNT + LEFT_NUC_MOTIF_COUNT

    plot = ggplot(model_coef_matrix, aes(x=position_value, y=base, fill=log_10_pdel)) +
        geom_tile() +
        theme_cowplot(font_family = 'Arial') + 
        xlab('Position') +
        ylab ('Base') +
        theme(text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 20))+
        geom_vline(xintercept = LEFT_NUC_MOTIF_COUNT + 0.5, size = 3.5, color = 'black') +
        guides(fill = guide_colourbar(barheight = 14)) +
        scale_fill_distiller(palette = 'PuOr', name = 'log10(probability of deletion)', limits = limits) +
        annotate("text", x = 0.35, y = 0.25, label = "5\'", size = 8) +  
        annotate("text", x = motif_length + 0.65, y = 0.25, label = "3\'", size = 8) +  
        coord_cartesian(ylim = c(1, 4), clip = "off")
    
    if (with_values == TRUE){
        plot = plot +
            geom_text(data = model_coef_matrix, aes(x = position_value, y = base, label = round(log_10_pdel, 3)))
    }

    return(plot)
}

plot_model_evaluation_loss_paracoord <- function(eval_data, loss_bound = NULL, color_palette = NULL, expand_var = 4, productivity_facet = TRUE) {
    if (is.factor(eval_data$long_loss_type)){
        types = levels(eval_data$long_loss_type)
        last = types[length(types)]
    } else {
        types = unique(eval_data$long_loss_type)
        last = types[length(types)]
    }
    
    label_data = eval_data[long_loss_type == last]

    # create plot
    require(ggrepel)
    plot = ggplot(eval_data) +
        geom_point(aes(y = log_loss, x = long_loss_type, color = long_model_type), size = 14)+
        geom_line(aes(y = log_loss, x = long_loss_type, group = long_model_type, color = long_model_type), size = 9, alpha = 0.8)+
        geom_text_repel(data = label_data, aes(y = log_loss, x = long_loss_type, label = long_model_type, color = long_model_type), nudge_x = 0.2, fontface = "bold", size = 13, direction = 'y', hjust = 0, point.padding = 1, max.overlaps = Inf, lineheight = 0.8) +
        theme_cowplot(font_family = 'Arial') + 
        xlab(' ') +
        ylab('Log loss\n') +
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) +
        theme(legend.position = 'none', text = element_text(size = 44), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 44), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))  +
        scale_x_discrete(expand = expansion(add = c(0.2, expand_var)))
    if (productivity_facet == TRUE){
        plot = plot + facet_wrap(~validation_productivity, ncol = 1)
    }

    if (!is.null(loss_bound)){
        plot = plot +
            ylim(loss_bound)
    }

    if (!is.null(color_palette)){
        plot = plot + scale_color_manual(values = color_palette)
    } else {
        plot = plot + scale_color_brewer(palette = 'Dark2')
    }

    return(plot)
}
