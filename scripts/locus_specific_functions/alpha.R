get_joining_gene_type <- function(){
    pairs = c('v_gene' = 'j_gene', 'j_gene' = 'v_gene', 'v-j_gene' = NULL)
    trims = c('v_trim' = 'j_trim', 'j_trim' = 'v_trim', 'v-j_trim' = NULL)
    inserts = c('v_gene' = 'vj_insert', 'j_gene' = 'vj_insert', 'v-j_gene' = 'vj_insert')
    return(list(gene_pairs = pairs, trim_pairs = trims, insert_pairs = inserts))
} 

pairs = get_joining_gene_type()

JOINING_GENE <<- pairs$gene_pairs[GENE_NAME]
JOINING_TRIM <<- pairs$trim_pairs[TRIM_TYPE]
JOINING_INSERT <<- pairs$insert_pairs[GENE_NAME]
CHAIN_TYPE <<- 'TRA'
