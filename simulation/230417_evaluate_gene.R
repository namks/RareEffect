# need to setwd first!
# setwd(path to gene_effect_all.txt)

library(data.table)
d <- fread("gene_effect_all.txt")

methods <- c("bayes", "bayes2", "linreg", "marginal", "ridge")

for (i in 1:2) {
    for (j in 1:2) {
        for (k in 1:2) {
            for (m in methods) {
                d2 <- d[which((d$method == m) & (d$ii == i) & (d$jj == j) & (d$kk == k)),]
                print(paste0(m, " scenario (", i, ", ", j, ", ", k, ") ", round(mean(d2$gene_effect), 4)))
            }
        }
    }
}

d3 <- d[which((d$method == "bayes2") & (d$ii == 2) & (d$jj == 2) & (d$kk == 2)),]
boxplot(gene_effect ~ gene, data = d3)
