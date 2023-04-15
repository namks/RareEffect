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
                    
                    # Bayesian
                    post_betas <- fread(paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/step2/single/bayes/sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt"), header=F)
                    lof <- betas_lof[betas_lof$V1 %in% post_betas$V1,]
                    lof2 <- left_join(lof, post_betas, by="V1")
                    mis <- betas_mis[betas_mis$V1 %in% post_betas$V1,]
                    mis2 <- left_join(mis, post_betas, by="V1")
                    syn <- betas_syn[betas_syn$V1 %in% post_betas$V1,]
                    syn2 <- left_join(syn, post_betas, by="V1")
                    all <- betas_all[betas_all$V1 %in% post_betas$V1,]
                    all2 <- left_join(all, post_betas, by="V1")

                    # Bayesian 2
                    post_betas_bayes2 <- fread(paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/step2/single/bayes2/sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt"), header=F)
                    lof3 <- left_join(lof2, post_betas_bayes2, by="V1")
                    mis3 <- left_join(mis2, post_betas_bayes2, by="V1")
                    syn3 <- left_join(syn2, post_betas_bayes2, by="V1")
                    all3 <- left_join(all2, post_betas_bayes2, by="V1")

                    # Linear regression
                    post_betas_linreg <- fread(paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/step2/single/linreg/sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt"), header=F)
                    lof4 <- left_join(lof3, post_betas_linreg, by="V1")
                    mis4 <- left_join(mis3, post_betas_linreg, by="V1")
                    syn4 <- left_join(syn3, post_betas_linreg, by="V1")
                    all4 <- left_join(all3, post_betas_linreg, by="V1")

                    # Marginal effect (linear regression)
                    post_betas_marginal <- fread(paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/step2/single/marginal/sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt"), header=F)
                    lof5 <- left_join(lof4, post_betas_marginal, by="V1")
                    mis5 <- left_join(mis4, post_betas_marginal, by="V1")
                    syn5 <- left_join(syn4, post_betas_marginal, by="V1")
                    all5 <- left_join(all4, post_betas_marginal, by="V1")

                    # Ridge regression
                    post_betas_ridge <- fread(paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/step2/single/ridge/sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt"), header=F)
                    lof6 <- left_join(lof5, post_betas_ridge, by="V1")
                    mis6 <- left_join(mis5, post_betas_ridge, by="V1")
                    syn6 <- left_join(syn5, post_betas_ridge, by="V1")
                    all6 <- left_join(all5, post_betas_ridge, by="V1")

                    # Write estimated betas
                    betas_lof_name <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/betas/betas_lof_sim", n, "_gene_", i, "_", ii, "_", jj, "_", kk, ".txt")
                    betas_mis_name <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/betas/betas_mis_sim", n, "_gene_", i, "_", ii, "_", jj, "_", kk, ".txt")
                    betas_syn_name <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/betas/betas_syn_sim", n, "_gene_", i, "_", ii, "_", jj, "_", kk, ".txt")
                    betas_all_name <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/betas/betas_all_sim", n, "_gene_", i, "_", ii, "_", jj, "_", kk, ".txt")

                    colnames(lof6) <- c("SNP", "True", "Bayes", "Bayes2", "Linreg", "Marginal", "Ridge")
                    colnames(mis6) <- c("SNP", "True", "Bayes", "Bayes2", "Linreg", "Marginal", "Ridge")
                    colnames(syn6) <- c("SNP", "True", "Bayes", "Bayes2", "Linreg", "Marginal", "Ridge")
                    colnames(all6) <- c("SNP", "True", "Bayes", "Bayes2", "Linreg", "Marginal", "Ridge")
                    write.table(lof6, betas_lof_name, row.names=F, col.names=F, quote=F)
                    write.table(mis6, betas_mis_name, row.names=F, col.names=F, quote=F)
                    write.table(syn6, betas_syn_name, row.names=F, col.names=F, quote=F)
                    write.table(all6, betas_all_name, row.names=F, col.names=F, quote=F)
                }

                print(paste("Scenario", ii, jj, kk, "for gene", i, "done."))
            }
        }
    }
}
