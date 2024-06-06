library(data.table)
library(mclogit)

d = data.table()
garden_count = 25
colors = c('red', 'yellow', 'blue')
species = c('dahlia', 'daisy', 'rose', 'marigold')

for (c in colors){
    for (s in species){
        h = sample(seq(10), 1)
        w = sample(seq(6), 1)
        r = sample(c('A', 'B', 'C', 'D'), 1)
        for (g in seq(garden_count)){
            ct = sample(seq(100), 1)
            temp = data.table(garden = g, color = c, species = s, count = ct, height = h, width = w, rating = r)
            d = rbind(d, temp)
        }
    }
}

d[, prob := count/sum(count), by = .(color, garden)]

fwrite(d, '/home/mrussel2/microhomology/jax_scripts/comparison_test_R/long_test_data.tsv', sep = '\t')

d$rating = factor(d$rating)
contrasts(d$rating) = contr.sum(4)

formula = formula(paste0('cbind(count, interaction(garden, color)) ~ height + width + rating')) 
lmodel = mclogit(formula = formula, data = d)  
