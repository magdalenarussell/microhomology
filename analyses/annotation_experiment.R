library(data.table)
library(vroom)

soto_vdj = fread('/fh/fast/matsen_e/shared/tcr-gwas/soto_alpha/vdjserver/pooled_processed.tsv')

igor_files = list.files('/fh/fast/matsen_e/shared/tcr-gwas/soto_alpha/igor', full.names = TRUE) 
soto_igor = as.data.table(vroom(igor_files))

adaptive_files = list.files('/fh/fast/matsen_e/shared/tcr-gwas/thymus_alpha/adaptive', full.names = TRUE)
adaptive = as.data.table(vroom(adaptive_files)) 

#filter adaptive
source('scripts/processing_functions/data_type_functions/adaptive.R')
LOCUS <<- 'TRA'
adaptive = convert_adaptive_style_to_imgt(adaptive)

# add type column
adaptive$type = 'adaptive'
soto_igor$type = 'soto igor'
soto_vdj$type = 'soto vdjserver'

cols = c('v_gene', 'j_gene', 'v_trim', 'j_trim', 'vj_insert', 'type')

# subset data
adaptive_subset = adaptive[v_trim >= 0 & v_trim <= 20 & j_trim >= 0 & j_trim <= 20 & !is.na(j_gene) & v_gene %like% 'TRA' & j_gene %like% 'TRA'][productive == 'nonproductive']
soto_igor_subset = soto_igor[v_trim >= 0 & v_trim <= 20 & j_trim >= 0 & j_trim <= 20 & !is.na(j_gene) & v_gene %like% 'TRA' & j_gene %like% 'TRA'] 
soto_vdj_subset = soto_vdj[v_trim >= 0 & v_trim <= 20 & j_trim >= 0 & j_trim <= 20 & !is.na(j_gene) & v_gene %like% 'TRA' & j_gene %like% 'TRA'][productive == 'nonproductive']

total = nrow(adaptive_subset)

# sample from soto data
allprob = data.table()
for (i in seq(100)){
    adaptive_sample = adaptive_subset[sample(.N,total,replace=TRUE)]
    soto_igor_sample = soto_igor_subset[sample(.N,total,replace=TRUE)]
    soto_vdj_sample = soto_vdj_subset[sample(.N,total,replace=TRUE)]

    together = rbind(adaptive_sample[, ..cols], soto_igor_sample[, ..cols], soto_vdj_sample[, ..cols])

    # fill in missing
    prob = together[, .N, by = .(v_gene, v_trim, type)]

    for (vg in unique(prob$v_gene)){
        for (vt in unique(prob$v_trim)){
            for (t in unique(prob$type)){
                temp = prob[v_gene == vg & v_trim == vt & type == t]
                if (nrow(temp) < 1){
                    new = data.table(v_gene = vg, v_trim = vt, type = t, N = 0)
                    prob= rbind(prob, new, fill = TRUE)
                }
            }
        }
    }

    prob[, total := sum(N), by = .(type, v_gene)]
    prob$iter = i
    allprob = rbind(prob, allprob)
}

allprob[, freq := N/total]
allprob[, all := sum(N), by = .(type, iter)]
allprob[, gene_freq := total/all]
allprob$v_gene = factor(allprob$v_gene, levels = unique(allprob[order(-gene_freq)]$v_gene))

allprob[, mean := mean(freq), by = .(v_gene, v_trim, type)]
allprob[, min := min(freq), by = .(v_gene, v_trim, type)] 
allprob[, max := max(freq), by = .(v_gene, v_trim, type)]

cols = c('v_gene', 'v_trim', 'type', 'mean', 'min', 'max')
subset = unique(allprob[, ..cols])

library(cowplot)
library(ggplot2)
plot = ggplot(subset) +
    geom_line(aes(x = v_trim, y = mean, group = type, color = type), size = 1.5)+
    geom_ribbon(aes(x = v_trim, ymin = min, ymax = max, fill = type), alpha = 0.4)+
    facet_wrap(~v_gene, ncol = 6)+
    xlab('V-gene trimming length') +
    ylab('Empirical probability') +
    theme_cowplot()+
    theme(text = element_text(size = 25), axis.text.x=element_text(size = 20), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + 
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) 

filepath = 'plots/annotation_experiment.pdf' 
ggsave(filepath, plot = plot, width = 30, height = 40, units = 'in', dpi = 750, device = cairo_pdf)
