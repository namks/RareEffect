library(data.table)
setwd("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/step2/gene/")

methods <- c("bayes", "bayes2", "linreg", "marginal", "ridge")

out <- NULL
for (m in methods) {
    for (ii in 1:2) {
        for (jj in 1:2) {
            for (kk in 1:2) {
                for (i in 1:10) {
                    for (n in 1:100) {
                        fname <- paste0(m, "/sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt")
                        if (file.exists(fname)) {
                            t <- fread(fname)
                            d <- c(m, n, i, ii, jj , kk, as.numeric(t[1, 1]))
                            out <- rbind(out, d)
                        }
                    }
                    # print(paste0("Gene-level effect size of gene ", i, " estimated by method ", m, " under scenario (", ii, ", ", jj, ", ", kk, ")"))
                    # print(paste0("Mean : ", round(mean(d, na.rm=T), 5)))
                    # print(paste0("Std. dev. : ", round(sd(d, na.rm=T), 5)))
                    print(paste(m, i))
                }
            }
        }
    }
}

colnames(out) <- c("method", "sim", "gene", "ii", "jj", "kk", "gene_effect")
write.table(out, "gene_effect_all.txt", row.names=F, quote=F)
