get_joining_gene_type <- function(){
    pairs = c('v_gene' = 'd_gene', 'j_gene' = 'd_gene')
    trims = c('v_trim' = 'd0_trim', 'j_trim' = 'd1_trim')
    inserts = c('v_gene' = 'vd_insert', 'j_gene' = 'dj_insert')
    return(list(gene_pairs = pairs, trim_pairs = trims, insert_pairs = inserts))
} 

pairs = get_joining_gene_type()

JOINING_GENE <<- pairs$gene_pairs[GENE_NAME]
JOINING_TRIM <<- pairs$trim_pairs[TRIM_TYPE]
JOINING_INSERT <<- pairs$insert_pairs[GENE_NAME]
CHAIN_TYPE <<- 'TRB'
