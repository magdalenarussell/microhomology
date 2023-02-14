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
sum_result = data.table()

index_list = c()
for (gene_name in unique(top_genes)){
    index = which(top_genes == gene_name)
    empirical_data = condensed_rep_data[get(GENE_NAME) == gene_name]
    sum_diverged = calculate_sum_of_abs_diffs(empirical_data, empirical_data)
    mean_sum_divergence = mean(sum_diverged$sum_abs_diff)

    subset = unique(empirical_data[[JOINING_GENE]])
    cols = c(JOINING_GENE, 'avg_paired_freq')
    subset_freq = unique(empirical_data[, ..cols])$avg_paird_freq

    # get gene sequence
    gene_seq = get_gene_sequence(whole_nucseqs, gene_name, gene_seq_length = 27, pnuc_count = 2)
    gene_seq_with_positions = get_plot_positions_for_gene_sequence(gene_seq)
        
    seq_text = 5
    joining_genes = whole_nucseqs[gene %in% subset]
    if (nrow(joining_genes) > 1){
        dists = get_similarity_distance(nt_count = 10, gene_type = toupper(substring(JOINING_GENE, 1, 1)), joining_genes)

        sum_diverged = merge(sum_diverged, dists, by.x = JOINING_GENE, by.y = 'gene', all.x = TRUE)

        clow = min(clow, min(empirical_data$dist - 0.01, na.rm = TRUE))
        chigh = max(chigh, max(empirical_data$dist + 0.01, na.rm = TRUE))
        index_list = c(index_list, index)
    } else {
        next
    }

    set.seed(55)
    boot = run_diff_bootstrap(subset, empirical_data, boot_count = 1000)
    sum_boot = boot$sum
    cols = c(GENE_NAME, 'bootstrap')
    mean_sum_divergence_boot = sum_boot[, mean(sum_abs_diff), by = cols]

    setnames(mean_sum_divergence_boot, 'V1', 'mean_divergence')

    sum_p = run_t_test(mean_sum_divergence_boot, mean_sum_divergence, nrow(sum_diverged)) 

    sum_result = rbind(sum_result, sum_p)

    temp_plot = ggplot() +
        geom_line(data = empirical_data, aes(x = get(TRIM_TYPE), y = p_trim_given_pair, group = get(JOINING_GENE)), size = 0.6, alpha = 0.5) +
        geom_vline(xintercept = 0, color = 'black', size = 2) +
        geom_text(data = sum_p, x = 20, y = 0.15, aes(label = paste0('p = ', signif(p, 4))), size = 9) +
        ggtitle(gene_name) +
        xlab('Number of trimmed nucleotides') +
        ylab('Probability') +
        theme_cowplot(font_family = 'Arial') + 
    theme(text = element_text(size = 25), axis.text.x=element_text(size = 20), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5), legend.position = 'bottom', legend.direction = 'horizontal', legend.justification="center") + 
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5)
 
    legend = get_legend(temp_plot)
    temp_plot = temp_plot + 
        xlab('') +
        ylab('') +
        theme(legend.position = "none", text = element_text(size = 20))

    assign(paste0('gene', index), temp_plot)
}

# combine all distributions
all = plot_grid(plotlist=mget(paste0("gene", index_list)), ncol = 5) 

# save plot
file_name = paste0(path, '/tra_', TRIM_TYPE, '_gallery_by_joining_gene_signif.pdf')
ggsave(file_name, plot = all, width = 35, height = 30, units = 'in', dpi = 750, device = cairo_pdf, limitsize = FALSE)

gene15 = gene15 +
    xlab('Number of trimmed nucleotides') +
    ylab('Probability') +
    theme(legend.position = "none", text = element_text(size = 18))

file_name = paste0(path, '/TRAV13-1*01_by_joining_gene_signif.pdf')
ggsave(file_name, plot = gene15, width = 9, height = 6, units = 'in', dpi = 750, device = cairo_pdf, limitsize = FALSE)


