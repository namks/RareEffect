# Install packages

if (!require(seqminer)) install.packages(seqminer)
if (!require(data.table)) install.packages(data.table)

# Load libraries

setwd("/media/leelabsg-storage0/kisung/RVPRS/simulation/script")
library(seqminer)
library(data.table)
source("RVPRS_function.R")

# Generate covariates and betas
set.seed(0)
n_sim <- 100

for (n in 1:n_sim) {
    groupfile_name <- "/media/leelabsg-storage0/kisung/dnanexus/group_files/UKBexome_all_chr.txt"
    gene_list <- c("ANGPTL3", "ANGPTL8", "APOA1", "CETP", "DOCK6", "FAM86C1", "LPL", "PLA2G12A", "PLIN1", "SMARCA4")

    for (i in seq_len(length(gene_list))) {
        gene_name <- as.character(gene_list[i])

        rare_plink_prefix <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/rare_variants/UKBB_200k_WES_", gene_name, "_WB")
        rare_plink_matrix <- load_plink(rare_plink_prefix)
        rare_plink_matrix_converted <- convert_missing(rare_plink_matrix, ref_first = FALSE)
        n_sample <- nrow(rare_plink_matrix_converted)

        # mat_by_func_anno = [G_lof, G_mis, G_syn]
        var_by_func_anno <- read_groupfile(groupfile_name, gene_name)
        mat_by_func_anno <- split_plink_matrix(rare_plink_matrix_converted,
            var_by_func_anno[[1]],
            var_by_func_anno[[2]],
            var_by_func_anno[[3]],
            mac_threshold = 1
        )

        # Genotype matrix for rare_variants
        lof_mat <- mat_by_func_anno[[1]][, c(1:(ncol(mat_by_func_anno[[1]]) - 1))]
        mis_mat <- mat_by_func_anno[[2]][, c(1:(ncol(mat_by_func_anno[[2]]) - 1))]
        syn_mat <- mat_by_func_anno[[3]][, c(1:(ncol(mat_by_func_anno[[3]]) - 1))]

        # Proportion of causal variants
        ## Scenario : (LoF, mis, syn) = (20%, 10%, 2%), (30%, 10%, 2%)
        lof_causal_rates <- c(0.2, 0.3)
        mis_causal_rates <- c(0.1, 0.1)
        syn_causal_rates <- c(0.02, 0.02)

        # Absolute effect sizes
        ## Scenarios : (LoF, mis, syn) = (0.3 log MAF, 0.15 log MAF, 0.15 log MAF), (0.5 log MAF, 0.25 log MAF, 0.25 log MAF)
        lof_effect_sizes <- c(0.3, 0.5)
        mis_effect_sizes <- c(0.15, 0.25)
        syn_effect_sizes <- c(0.15, 0.25)

        # Effect directions
        ## Scenarios : (LoF, mis, syn) = (100%, 80%, 50%), (100%, 100%, 100%)
        lof_direction <- c(1, 1)
        mis_direction <- c(0.8, 1)
        syn_direction <- c(0.5, 1)

        # MAC by functional annotation
        lof_MAC <- colSums(lof_mat)
        mis_MAC <- colSums(mis_mat)
        syn_MAC <- colSums(syn_mat)

        # Generate causal vector
        betas_lof <- rep(0, ncol(lof_mat))
        betas_mis <- rep(0, ncol(mis_mat))
        betas_syn <- rep(0, ncol(syn_mat))

        for (ii in 1:2) {
            for (jj in 1:2) {
                for (kk in 1:2) {
                    # Generate beta for LoF
                    for (j in seq_len(ncol(lof_mat))) {
                        MAF <- lof_MAC[j] * 2 / n_sample
                        if (lof_MAC[j] <= 10) {
                            betas_lof[j] <- rbinom(1, 1, 3 * lof_causal_rates[ii]) * abs(lof_effect_sizes[jj] * log(MAF, base = 10))
                        } else {
                            betas_lof[j] <- rbinom(1, 1, lof_causal_rates[ii]) * abs(lof_effect_sizes[jj] * log(MAF, base = 10))
                        }
                    }

                    # Generate beta for mis
                    for (j in seq_len(ncol(mis_mat))) {
                        MAF <- mis_MAC[j] * 2 / n_sample
                        if (mis_MAC[j] <= 10) {
                            betas_mis[j] <- rbinom(1, 1, 3 * mis_causal_rates[ii]) * abs(mis_effect_sizes[jj] * log(MAF, base = 10)) * (2 * (rbinom(1, 1, mis_direction[kk]) - 0.5))
                        } else {
                            betas_mis[j] <- rbinom(1, 1, mis_causal_rates[ii]) * abs(mis_effect_sizes[jj] * log(MAF, base = 10)) * (2 * (rbinom(1, 1, mis_direction[kk]) - 0.5))
                        }
                    }

                    # Generate beta for syn
                    for (j in seq_len(ncol(syn_mat))) {
                        MAF <- syn_MAC[j] * 2 / n_sample
                        if (syn_MAC[j] <= 10) {
                            betas_syn[j] <- rbinom(1, 1, 3 * syn_causal_rates[ii]) * abs(syn_effect_sizes[jj] * log(MAF, base = 10)) * (2 * (rbinom(1, 1, syn_direction[kk]) - 0.5))
                        } else {
                            betas_syn[j] <- rbinom(1, 1, syn_causal_rates[ii]) * abs(syn_effect_sizes[jj] * log(MAF, base = 10)) * (2 * (rbinom(1, 1, syn_direction[kk]) - 0.5))
                        }
                    }

                    # Phenotype generation (rare-variant part)
                    pheno_lof <- lof_mat %*% betas_lof
                    pheno_mis <- mis_mat %*% betas_mis
                    pheno_syn <- syn_mat %*% betas_syn

                    pheno_rv <- pheno_lof + pheno_mis + pheno_syn

                    betas_lof_name <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/data/betas/betas_lof_sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt")
                    betas_mis_name <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/data/betas/betas_mis_sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt")
                    betas_syn_name <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/data/betas/betas_syn_sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt")

                    out <- cbind(rownames(pheno_rv), pheno_rv)
                    colnames(out) <- c("IID", "pheno_rv")
                    
                    outname <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/data/pheno/pheno_sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt")

                    write.table(betas_lof, betas_lof_name, row.names=F, col.names=F, quote=F)
                    write.table(betas_mis, betas_mis_name, row.names=F, col.names=F, quote=F)
                    write.table(betas_syn, betas_syn_name, row.names=F, col.names=F, quote=F)
                    write.table(out, outname, row.names=F, col.names=T, quote=F)                    
                }
            }
        }
    }
    print(paste0("Simulation ", n, " completed"))
}
