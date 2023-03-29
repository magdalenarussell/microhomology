library(data.table)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

output_directory = args[1]
final_output_directory = args[2]

source('config/config.R')
source('scripts/igor_scripts/igor_processing_functions.R')

output = fread(file.path(output_directory, 'foo_output/parsed_sampled_output.csv'))
alignments = fread(file.path(output_directory, 'aligns/foo_indexed_CDR3s.csv'))

together = merge(alignments, output)

cols = c('seq_index', 'scenario_rank', 'scenario_proba_cond_seq', 'CDR3nt', 'CDR3aa', 'v_gene_call', 'j_gene_call', 'v_trim', 'j_trim', 'vj_insert', 'vj_insert_nucs')

together = together[, ..cols]

setnames(together, 'v_gene_call', 'v_gene')
setnames(together, 'j_gene_call', 'j_gene')
setnames(together, 'CDR3nt', 'cdr3_nucseq')
setnames(together, 'CDR3aa', 'cdr3')

# extract just the gene names
together[, v_gene := get_gene_from_long_igor_column(v_gene, sep = "\\;")]
together[, j_gene := get_gene_from_long_igor_column(j_gene, sep = "\\;")]

#assign productivity status
together[cdr3 == '', productive := FALSE]
together[cdr3 != '', productive := TRUE]

#assign pnucs to their own category
for (gene in c('v', 'j')){
    trim = paste0(gene, '_trim')
    pnuc = paste0(gene, '_pnuc')
    together[, paste(pnuc) := 0]
    together[get(trim) < 0, paste(pnuc) := -1*(get(trim))]
    together[get(trim) < 0, paste(trim) := 0]
}

file_name = create_final_output_location(output_directory, final_output_directory)
fwrite(together, file_name, sep = '\t')
