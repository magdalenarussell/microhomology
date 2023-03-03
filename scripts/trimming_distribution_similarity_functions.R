pairwise_diffs <- function(model_predictions){
    # pairwise sum of abs diffs
    new_dt = data.table() 
    for (v in unique(model_predictions[[GENE_NAME]])){
        temp = model_predictions[get(GENE_NAME) == v]
        for (j1 in unique(temp[[JOINING_GENE]])){
            for (j2 in unique(temp[[JOINING_GENE]])){
                if (j1 == j2){
                    next
                }
                j1temp = temp[get(JOINING_GENE) == j1]
                j2temp = temp[get(JOINING_GENE) == j2]
                tog = merge(j1temp, j2temp, by = c('v_trim', 'v_gene'))
                new_dt = rbind(new_dt, tog)
            }

        }
    }
 
    new_dt[, abs_diff := abs(avg_prob.x - avg_prob.y)]
    pairwise = new_dt[, sum(abs_diff), by = .(v_gene, j_gene.x, j_gene.y)]
    setnames(pairwise, 'V1', 'sum_abs_diff')
    return(pairwise)
}

cluster_diffs <- function(model_predictions, cluster_count = 3){
    new_dt = data.table()
    for (v in unique(model_predictions[[GENE_NAME]])){
        temp = model_predictions[get(GENE_NAME) == v]
        wider = temp %>% 
            pivot_wider(names_from = v_trim, values_from = avg_prob) %>%
            as.data.table()
        km = kmeans(wider[, -c('j_gene', 'v_gene')], centers = cluster_count, nstart = 25)
        wider = cbind(wider, km$cluster) 
        setnames(wider, 'V2', 'cluster')

        temp = merge(temp, wider[, c('v_gene', 'j_gene', 'cluster')])
        new_dt = rbind(new_dt, temp)
    }
    return(new_dt)
}
