setwd("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/betas/")

beta_all <- NULL
func_anno <- c("all", "lof", "mis", "syn")

for (i in 1:10) {
    for (j in 1:4) {
        for (ii in 1:2) {
            for (jj in 1:2) {
                for (kk in 1:2) {
                    for (n in 1:100) {
                        fname <- paste0("betas_", func_anno[j], "_sim", n, "_gene_", i, "_", ii, "_", jj, "_", kk, ".txt")
                        try({ 
                            d <- read.table(fname)
                            if (nrow(d) > 0) {
                                temp <- cbind(func_anno[j], n, i, ii, jj, kk, d)
                                beta_all <- rbind(beta_all, temp)
                            }
                        })
                    }
                    print(paste0("Done for simulations with gene ", i, " and functional annotation ", j, " and ", ii, ", ", jj, ", ", kk))
                }
            }
        }
    }
}

colnames(beta_all) <- c("func_anno", "n", "gene", "ii", "jj", "kk", "variant", "true_beta", "bayes", "bayes2", "linear", "marginal", "ridge")
outname <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/betas/beta_all_v2.txt")
write.table(beta_all, outname, row.names=F, quote=F)
