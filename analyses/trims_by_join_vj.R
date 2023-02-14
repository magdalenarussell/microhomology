source('config/config.R')

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(cowplot)

TRIM_TYPE <<- 'v_trim' 
JOINING_GENE <<- 'j_gene'
PRODUCTIVITY <<- 'nonproductive' 
NCPU <<- as.numeric(10)
GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')

path = paste0(PROJECT_PATH, '/plots/trim_by_join')
dir.create(path, recursive = TRUE)

source(paste0(PROJECT_PATH, '/scripts/data_processing_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/plotting_functions.R'))

# get whole germline sequences and their distances from one another
whole_nucseqs = get_whole_nucseqs()
dists = get_similarity_distance(nt_count = 10, gene_type = toupper(substring(JOINING_GENE, 1, 1)), whole_nucseqs)

# Compile repertoire data for all subjects
rep_data = read_all_data(directory = TRA_REPERTOIRE_DATA)
rep_data = convert_adaptive_style_to_imgt(rep_data) 
rep_data_subset = filter_data(rep_data)
condensed_rep_data = condense_data(rep_data_subset) 

# order genes by frequency
top_genes = unique(condensed_rep_data[order(-avg_paired_freq)][[GENE_NAME]])

# Create a trimming distribution plot for each gene
all_plots = c()

max_prob = 1
count = 0
clow = 1
chigh = 0
index_list = c()
for (gene_name in unique(top_genes)){
    index = which(top_genes == gene_name)
    empirical_data = condensed_rep_data[get(GENE_NAME) == gene_name]
    diverged = calculate_sum_of_abs_diffs(empirical_data)

    subset = unique(empirical_data[[JOINING_GENE]])
    # get gene sequence
    gene_seq = get_gene_sequence(whole_nucseqs, gene_name, gene_seq_length = 27, pnuc_count = 2)
    gene_seq_with_positions = get_plot_positions_for_gene_sequence(gene_seq)
        
    seq_text = 5
    joining_genes = whole_nucseqs[gene %in% subset]
    if (nrow(joining_genes) > 1){
        dists = get_similarity_distance(nt_count = 10, gene_type = toupper(substring(JOINING_GENE, 1, 1)), joining_genes)

        empirical_data = merge(empirical_data, dists, by.x = JOINING_GENE, by.y = 'gene', all.x = TRUE)
        diverged = merge(diverged, dists, by.x = JOINING_GENE, by.y = 'gene', all.x = TRUE)

        count = count + 1
        clow = min(clow, min(empirical_data$dist - 0.01, na.rm = TRUE))
        chigh = max(chigh, max(empirical_data$dist + 0.01, na.rm = TRUE))
        index_list = c(index_list, index)
    } else {
        next
    }

    empirical_data = empirical_data[order(dist)]
    empirical_data[[JOINING_GENE]] = factor(empirical_data[[JOINING_GENE]], levels = unique(empirical_data[[JOINING_GENE]]))

    temp_plot = ggplot() +
        geom_line(data = empirical_data, aes(x = get(TRIM_TYPE), y = p_trim_given_pair, group = get(JOINING_GENE), color = dist), size = 0.6) +
        geom_vline(xintercept = 0, color = 'black', size = 2) +
        geom_text(data = gene_seq_with_positions, y = max_prob, aes(x = position, label = base), size = seq_text, color = 'black') +
        geom_text(y = max_prob, aes(x = -2.1), label = '3\'- ', size = seq_text - 1) +
        geom_text(y = max_prob, aes(x = 24 + 0.1), label = ' -5\'', size = seq_text - 1) +
        ggtitle(gene_name) +
        xlab('Number of trimmed nucleotides') +
        ylab('Probability') +
        theme_cowplot(font_family = 'Arial') + 
    theme(text = element_text(size = 25), axis.text.x=element_text(size = 20), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5), legend.position = 'bottom', legend.direction = 'horizontal', legend.justification="center") + 
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5)+
        scale_color_gradient(limits = c(clow, chigh),low = "yellow", high = "red", na.value = 'gray')+
        guides(colour = guide_colorbar(barwidth = 35))+
        labs(color = 'Avg pairwise hamming\ndistance of first 10 nt')
 
    legend = get_legend(temp_plot)
    temp_plot = temp_plot + 
        xlab('') +
        ylab('') +
        theme(legend.position = "none", text = element_text(size = 20), plot.margin = margin(l = -0.1, b = -0.3, unit = 'cm'))

    assign(paste0('gene', index), temp_plot)

    # plot divergence versus difference

    rel = lm(dist ~ sum_abs_diff, data = diverged)
    relation = data.table(r2 = round(summary(rel)$r.squared, 3), slope = round(coef(rel)[['sum_abs_diff']], 3), intercept = round(coef(rel)[['(Intercept)']], 3))

    temp_plot2 = ggplot(diverged) +
        geom_point(aes(x = sum_abs_diff, y = dist), alpha = 0.6, size = 5) +
        ggtitle(gene_name) +
        geom_text(data = relation, x = 0.08, y = 0.4, aes(label = paste0('R^2 = ', r2)), size = 9) +
        xlab('Sum of absolute differences\nfrom mean distribution')+
        ylab('Average pairwsie hamming\ndistance of the first 10 nt')+
        theme_cowplot(font_family = 'Arial') + 
    theme(text = element_text(size = 25), axis.text.x=element_text(size = 20), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + 
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) +
        xlim(0.04, 0.23) +
        ylim(0.35, 0.78)

    assign(paste0('div', index), temp_plot2)
}

# combine all distributions
all = plot_grid(plotlist=mget(paste0("gene", index_list)), ncol = 5) 
all = plot_grid(all, NULL, legend, NULL, ncol = 1, rel_heights = c(1, 0.05, 0.1, 0.05))

all2 = plot_grid(plotlist=mget(paste0("div", index_list)), ncol = 5) 

# save plot
file_name = paste0(path, '/tra_', TRIM_TYPE, '_gallery_by_joining_gene.pdf')
ggsave(file_name, plot = all, width = 35, height = 30, units = 'in', dpi = 750, device = cairo_pdf, limitsize = FALSE)

file_name2 = paste0(path, '/tra_', TRIM_TYPE, '_divergence_by_joining_gene.pdf')
ggsave(file_name2, plot = all2, width = 35, height = 46, units = 'in', dpi = 750, device = cairo_pdf, limitsize = FALSE)

gene15 = gene15 +
 scale_color_gradient(limits = c(0.475, 0.74),low = "yellow", high = "red", na.value = 'gray') +
 theme(text = element_text(size = 18), axis.text.x=element_text(size = 15), axis.text.y = element_text(size = 15), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5), legend.position = 'bottom', legend.direction = 'horizontal', legend.justification="center") +
guides(colour = guide_colorbar(barwidth = 20))

legend1 = get_legend(gene15)

gene15 = gene15 +
    xlab('Number of trimmed nucleotides') +
    ylab('Probability') +
    theme(legend.position = "none", text = element_text(size = 18))


gene15p = plot_grid(gene15, NULL, legend1, NULL, ncol = 1, rel_heights = c(1, 0.05, 0.1, 0.05))
file_name = paste0(path, '/TRAV13-1*01_by_joining_gene.pdf')
ggsave(file_name, plot = gene15p, width = 9, height = 6, units = 'in', dpi = 750, device = cairo_pdf, limitsize = FALSE)

div15 = div15 +
        xlab('Sum of absolute differences\nfrom mean distribution')+
        ylab('Average pairwsie hamming\ndistance of the first 10 nt') +
        theme(legend.position = "none", text = element_text(size = 18))

file_name = paste0(path, '/TRAV13-1*01_by_joining_gene_divergence.pdf')
ggsave(file_name, plot = div15, width = 9, height = 8, units = 'in', dpi = 750, device = cairo_pdf, limitsize = FALSE)

