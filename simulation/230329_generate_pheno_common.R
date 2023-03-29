# Install packages

if (!require(seqminer)) install.packages(seqminer)
if (!require(data.table)) install.packages(data.table)

# Load libraries

setwd("/media/leelabsg-storage0/kisung/RVPRS/simulation/script")
library(seqminer)
library(data.table)
source("RVPRS_function.R")

plink_prefix <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/common_variants/common_variants_WB")
plink_matrix <- load_plink(plink_prefix)
plink_matrix_converted <- convert_missing(plink_matrix, ref_first = FALSE)
print("PLINK file for common variants loaded")

# Generate covariates and betas
set.seed(0)
n_sim <- 100
n_sample <- nrow(plink_matrix_converted)

for (n in 1:n_sim) {
    n_common <- ncol(plink_matrix_converted)
    cov1 <- rbinom(n_sample, 1, 0.1)
    cov2 <- rnorm(n_sample, 0, 1)
    betas <- rnorm(n_common, 0, 1 / n_common)
    noise <- rnorm(n_sample, 0, 1)

    # Phenotype generation (Common variant part)
    pheno_common <- plink_matrix_converted %*% betas

    out <- cbind(rownames(pheno_common), pheno_common, cov1, cov2, noise)
    colnames(out)[1] <- "IID"
    colnames(out)[2] <- "pheno_common"
    outname <- paste0("/media/leelabsg-storage0/kisung/RVPRS/simulation/data/pheno/pheno_sim", n, ".txt")

    write.table(out, outname, row.names=F, col.names=T, quote=F)
    print(paste("Simulation", n, "completed"))
}
