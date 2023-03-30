setwd("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/step1")
library(sparsesvd)
library(dplyr)
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
                
                # Run analysis be gene
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

                    system.time({
                        post_beta <- fast_lmm(G = G_reordered, Y = y_tilde[,2])
                    })
                    lof_ncol <- ncol(mat_by_func_anno[[1]])
                    gene_effect_size <- calc_gene_effect_size(G = G_reordered, lof_ncol = lof_ncol, post_beta = post_beta, Y = y_tilde[,2])

                    single_effect <- cbind(colnames(G), post_beta[,1])
                    
                    single_outname <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/step2/single/sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt")
                    gene_outname <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/step2/gene/sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt")

                    write.table(single_effect, single_outname, row.names=F, col.names=F, quote=F)
                    write.table(gene_effect_size, gene_outname, row.names=F, col.names=F, quote=F)
                }
            }
        }
    }
    print(paste0("Simulation ", n, " completed."))
}
