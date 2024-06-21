get_joining_gene_type <- function(){
    pairs = c('d_gene' = 'j_gene', 'j_gene' = 'd_gene', 'd-j_gene' = NULL)
    trims = c('d3_trim' = 'j_trim', 'j_trim' = 'd3_trim', 'd3-j_trim' = NULL)
    inserts = c('d_gene' = 'dj_insert', 'j_gene' = 'dj_insert', 'd-j_gene' = 'dj_insert')
    return(list(gene_pairs = pairs, trim_pairs = trims, insert_pairs = inserts))
} 


pairs = get_joining_gene_type()

JOINING_GENE <<- pairs$gene_pairs[GENE_NAME]
JOINING_TRIM <<- pairs$trim_pairs[TRIM_TYPE]
JOINING_INSERT <<- pairs$insert_pairs[GENE_NAME]
CHAIN_TYPE <<- 'IGH'
CHAIN_SUBTYPE <<- 'B'
JUNCTION_TYPE <<- 'VDJ'
SUB_JUNCTION_TYPE <<- 'DJ'
