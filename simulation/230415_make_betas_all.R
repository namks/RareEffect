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
                        d <- read.table(fname)
                        temp <- cbind(func_anno[j], n, i, ii, jj, kk, d)
                        beta_all <- rbind(beta_all, temp)
                    }
                }
            }
        }
    }
}

colnames(beta_all) <- c("func_anno", "n", "gene", "ii", "jj", "kk", "variant", "true_beta", "bayes", "bayes2", "linear", "marginal", "ridge")
outname <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/betas/beta_all_v2.txt")
write.table(beta_all, outname, row.names=F, quote=F)
