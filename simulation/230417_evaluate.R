setwd("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/betas/")
library(data.table)

d <- fread("beta_all_v2.txt")
d[is.na(d)] <- 0
func_anno <- c("all", "lof", "mis", "syn")
methods <- c("bayes", "bayes2", "linear", "marginal", "ridge")

rmse <- function(x, y) {
    return(sqrt(mean((x - y) ^ 2)))
}

# evaluation
print("Aggregated by all genes and all scenario...")
for (f in func_anno) {
    d2 <- d[which(d$func_anno == f),]
    print(paste0("functional annotation: ", f))
    for (m in methods) {
        print(paste0("Method ", m, " RMSE = ", round(rmse(d2$true_beta, unlist(d2[, ..m])), 4)))
    }
}

print("Evaluated by scenario")
for (i in 1:2) {
    for (j in 1:2) {
        for (k in 1:2) {
            for (f in func_anno) {
                print(paste0("Scenario (", i, ", ", j, ", ", k, ")"))
                d2 <- subset(d, (d$func_anno == f) & (d$ii == i) & (d$jj == j) & (d$kk == k))
                print(paste0("functional annotation: ", f))
                for (m in methods) {
                    print(paste0("Method ", m, " RMSE = ", round(rmse(d2$true_beta, unlist(d2[, ..m])), 4)))
                }
            }
        }
    }
}

print("Evaluated by gene")
for (i in 1:10) {
    for (f in func_anno) {
        print(paste0("Gene ", i))
        d2 <- subset(d, (d$func_anno == f) & (d$gene == i))
        print(paste0("functional annotation: ", f))
        for (m in methods) {
            print(paste0("Method ", m, " RMSE = ", round(rmse(d2$true_beta, unlist(d2[, ..m])), 4)))
        }
    }
}
