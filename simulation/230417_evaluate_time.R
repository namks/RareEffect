setwd("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/step2/time")
library(data.table)

methods <- c("bayes", "bayes2", "linreg", "marginal", "ridge")
d <- NULL
for (m in methods) {
    for (n in 1:100) {
        for (i in 1:10) {
            for (ii in 1:2) {
                for (jj in 1:2) {
                    for (kk in 1:2) {
                        fname <- paste0(m, "/sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt")
                        if (file.exists(fname)) {
                            t <- fread(fname, header=T)
                            t2 <- c(m, n, i, ii, jj, kk, as.numeric(t[1, 3]))
                            d <- rbind(d, t2)
                        }
                    }
                }
            }
        }
        print(paste0("Simulation ", n, " completed"))
    }
}

colnames(d) <- c("method", "sim", "gene", "ii", "jj", "kk", "time")
write.table(d, "time_all.txt", row.names=F, quote=F)

d <- fread("time_all.txt")
# All gene
print(paste0("Avg. time for all gene"))
for (m in methods) {
    d2 <- d[which(d$method == m),]
    print(paste0("Method ", m, " average time : ", round(mean(d2$time, na.rm=T), 4)))
}

# For each gene
for (i in 1:10) {
    print(paste0("Avg. time for gene ", i))
    for (m in methods) {
        d2 <- d[which((d$method == m) & (d$gene == i)),]
        print(paste0("Method ", m, " average time : ", round(mean(d2$time, na.rm=T), 4)))
    }
}
