setwd("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/step1")
library(sparsesvd)
library(dplyr)
library(glmnet)
source("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/RVPRS_function.R")

n_sim <- 100
groupfile_name <- "/media/leelabsg-storage0/kisung/dnanexus/group_files/UKBexome_all_chr.txt"
gene_list <- c("ANGPTL3", "ANGPTL8", "APOA1", "CETP", "DOCK6", "FAM86C1", "LPL", "PLA2G12A", "PLIN1", "SMARCA4")

for (n in 1:n_sim) {
    for (ii in 1:2) {
        for (jj in 1:2) {
            for (kk in 1:2) {
                fname <- paste0("step1_sim", n, "_", ii, "_", jj, "_", kk, ".rda")
                load(fname)

                y_tilde <- cbind(as.numeric(modglmm$sampleID), modglmm$residuals)
                
                # Run analysis by gene
                for (i in 1:length(gene_list)) {
                    plink_prefix <- paste0("/media/leelabsg-storage0/kisung/RVPRS/plink/UKBB_200k_WES_", gene_list[i], "_WB")
                    plink_matrix <- load_plink(plink_prefix)
                    plink_matrix_converted <- convert_missing(plink_matrix, ref_first = FALSE)

                    var_by_func_anno <- read_groupfile(groupfile_name, gene_list[i])
                    mat_by_func_anno <- split_plink_matrix(plink_matrix_converted,
                        var_by_func_anno[[1]],
                        var_by_func_anno[[2]],
                        var_by_func_anno[[3]],
                        mac_threshold = 10
                    )

                    G <- cbind(mat_by_func_anno[[1]], mat_by_func_anno[[2]], mat_by_func_anno[[3]])
                    G_reordered <- G[match(y_tilde[,1], rownames(G)),]

                    # Linear regression
                    time_elapsed <- system.time({
                        # Linear regression
                        m1 <- lm(y_tilde[,2] ~ G_reordered)
                        post_beta <- m1$coefficients[2:length(m1$coefficients)]
                    })
                    lof_ncol <- ncol(mat_by_func_anno[[1]])
                    gene_effect_size <- calc_gene_effect_size(G = G_reordered, lof_ncol = lof_ncol, post_beta = post_beta, Y = y_tilde[,2])

                    single_effect <- cbind(colnames(G), post_beta)
                    
                    time_outname <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/step2/time/linreg/sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt")
                    single_outname <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/step2/single/linreg/sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt")
                    gene_outname <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/step2/gene/linreg/sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt")

                    write.table(t(time_elapsed), time_outname, row.names=F, col.names=T, quote=F)
                    write.table(single_effect, single_outname, row.names=F, col.names=F, quote=F)
                    write.table(gene_effect_size, gene_outname, row.names=F, col.names=F, quote=F)

                    # Marginal effect
                    post_beta <- rep(0, ncol(G_reordered))
                    time_elapsed <- system.time({
                        # Linear regression by each variant (take marginal effect)
                        for (j in 1:ncol(G_reordered)) {
                            m1 <- lm(y_tilde[,2] ~ G_reordered[,j])
                            post_beta[j] <- m1$coefficients[2]
                        }
                    })
                    lof_ncol <- ncol(mat_by_func_anno[[1]])
                    gene_effect_size <- calc_gene_effect_size(G = G_reordered, lof_ncol = lof_ncol, post_beta = post_beta, Y = y_tilde[,2])

                    single_effect <- cbind(colnames(G), post_beta)
                    
                    time_outname <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/step2/time/marginal/sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt")
                    single_outname <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/step2/single/marginal/sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt")
                    gene_outname <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/step2/gene/marginal/sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt")

                    write.table(t(time_elapsed), time_outname, row.names=F, col.names=T, quote=F)
                    write.table(single_effect, single_outname, row.names=F, col.names=F, quote=F)
                    write.table(gene_effect_size, gene_outname, row.names=F, col.names=F, quote=F)

                    # Ridge regression
                    time_elapsed <- system.time({
                        # Ridge regression
                        cv_fit <- cv.glmnet(G_reordered, y_tilde[,2], alpha = 0)
                        opt_lambda <- cv_fit$lambda.min

                        m1 <- glmnet(G_reordered, y_tilde[,2], alpha = 0, lambda = opt_lambda)
                        post_beta <- coef(m1)[2:nrow(coef(m1))]
                    })
                    lof_ncol <- ncol(mat_by_func_anno[[1]])
                    gene_effect_size <- calc_gene_effect_size(G = G_reordered, lof_ncol = lof_ncol, post_beta = post_beta, Y = y_tilde[,2])

                    single_effect <- cbind(colnames(G), post_beta)
                    
                    time_outname <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/step2/time/ridge/sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt")
                    single_outname <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/step2/single/ridge/sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt")
                    gene_outname <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/step2/gene/ridge/sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt")

                    write.table(t(time_elapsed), time_outname, row.names=F, col.names=T, quote=F)
                    write.table(single_effect, single_outname, row.names=F, col.names=F, quote=F)
                    write.table(gene_effect_size, gene_outname, row.names=F, col.names=F, quote=F)
                }
            }
        }
    }
    print(paste0("Simulation ", n, " completed."))
}
