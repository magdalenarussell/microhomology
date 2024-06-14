get_joining_gene_type <- function(){
    pairs = c('v_gene' = 'd_gene', 'd_gene' = 'j_gene', 'v-d_gene' = NULL)
    trims = c('d5_trim' = 'v_trim', 'v_trim' = 'd5_trim', 'v-d_trim' = NULL)
    inserts = c('d_gene' = 'vd_insert', 'v_gene' = 'vd_insert', 'v-d_gene' = 'vd_insert')
    return(list(gene_pairs = pairs, trim_pairs = trims, insert_pairs = inserts))
}

pairs = get_joining_gene_type()

JOINING_GENE <<- pairs$gene_pairs[GENE_NAME]
JOINING_TRIM <<- pairs$trim_pairs[TRIM_TYPE]
JOINING_INSERT <<- pairs$insert_pairs[GENE_NAME]
CHAIN_TYPE <<- 'IGH'
JUNCTION_TYPE <<- 'VDJ'

get_gene_order <- function(gene_type){
    return(c('v_gene', 'd_gene'))
}

get_trim_order <- function(trim_type){
    final = c('v_trim', 'd5_trim')
    if (trim_type %like% 'adjusted_mh'){
        final = paste0(final, '_adjusted_mh')
    }
    return(final)
}


