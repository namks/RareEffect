library(data.table)
library(dplyr)
n_sim <- 100
setwd("/media/leelabsg-storage0/kisung/RVPRS/simulation/data/betas")

for (i in 1:10) {
    lof_var <- fread(paste0("lof_gene", i, ".txt"), header=F)
    mis_var <- fread(paste0("mis_gene", i, ".txt"), header=F)
    syn_var <- fread(paste0("syn_gene", i, ".txt"), header=F)
    for (ii in 1:2) {
        for (jj in 1:2) {
            for (kk in 1:2) {
                cor_lof <- rep(0, n_sim)
                cor_mis <- rep(0, n_sim)
                cor_syn <- rep(0, n_sim)
                cor_all <- rep(0, n_sim)

                for (n in 1:n_sim) {
                    betas_lof <- fread(paste0("betas_lof_sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt"), header=F)
                    betas_mis <- fread(paste0("betas_mis_sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt"), header=F)
                    betas_syn <- fread(paste0("betas_syn_sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt"), header=F)
                    betas_lof <- cbind(lof_var, betas_lof)
                    betas_mis <- cbind(mis_var, betas_mis)
                    betas_syn <- cbind(syn_var, betas_syn)
                    colnames(betas_lof) <- c("V1", "V2")
                    colnames(betas_mis) <- c("V1", "V2")
                    colnames(betas_syn) <- c("V1", "V2")

                    betas_all <- rbind(betas_lof, betas_mis, betas_syn)

                    post_betas <- fread(paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/step2/single/sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt"), header=F)
                    lof <- betas_lof[betas_lof$V1 %in% post_betas$V1,]
                    lof2 <- left_join(lof, post_betas, by="V1")
                    mis <- betas_mis[betas_mis$V1 %in% post_betas$V1,]
                    mis2 <- left_join(mis, post_betas, by="V1")
                    syn <- betas_syn[betas_syn$V1 %in% post_betas$V1,]
                    syn2 <- left_join(syn, post_betas, by="V1")
                    all <- betas_all[betas_all$V1 %in% post_betas$V1,]
                    all2 <- left_join(all, post_betas, by="V1")
                    
                    cor_lof[n] <- cor(lof2$V2.x, lof2$V2.y)
                    cor_mis[n] <- cor(mis2$V2.x, mis2$V2.y)
                    cor_syn[n] <- cor(syn2$V2.x, syn2$V2.y)
                    cor_all[n] <- cor(all2$V2.x, all2$V2.y)
                }

                cor_lof_name <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/corr/corr_lof_", ii, "_", jj, "_", kk, ".txt")
                cor_mis_name <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/corr/corr_mis_", ii, "_", jj, "_", kk, ".txt")
                cor_syn_name <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/corr/corr_syn_", ii, "_", jj, "_", kk, ".txt")
                cor_all_name <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/corr/corr_all_", ii, "_", jj, "_", kk, ".txt")

                write.table(cor_lof, cor_lof_name, row.names=F, col.names=F, quote=F)
                write.table(cor_mis, cor_mis_name, row.names=F, col.names=F, quote=F)
                write.table(cor_syn, cor_syn_name, row.names=F, col.names=F, quote=F)
                write.table(cor_all, cor_all_name, row.names=F, col.names=F, quote=F)

                print(paste("Scenario", ii, jj, kk, "done."))
            }
        }
    }
}
