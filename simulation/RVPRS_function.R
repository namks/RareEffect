# Load packages

if (!require(seqminer)) install.packages(seqminer)
if (!require(data.table)) install.packages(data.table)
if (!require(tibble)) install.packages(tibble)
if (!require(dplyr)) install.packages(dplyr)
if (!require(Matrix)) install.packages(Matrix)
if (!require(sparsesvd)) install.packages(sparsesvd)

library(seqminer)
library(data.table)
library(tibble)
library(dplyr)
library(Matrix)
library(sparsesvd)

# Load PLINK file
load_plink <- function(plink_prefix) {
    plink_obj <- seqminer::openPlink(plink_prefix)
    marker_index <- seq(nrow(plink_obj$bim))
    sample_index <- seq(nrow(plink_obj$fam))
    plink_matrix <- seqminer::readPlinkToMatrixByIndex(plink_prefix, sample_index, marker_index)
    colnames(plink_matrix) <- plink_obj$bim$V2

    return(plink_matrix)
}

# Read Group file and split variants by functional annotations
read_groupfile <- function(groupfile_name, gene_name) {
    groupfile <- file(groupfile_name, "r")
    line <- 0
    var <- NULL
    anno <- NULL

    while (TRUE) {
        line <- line + 1
        marker_group_line <- readLines(groupfile, n = 1)

        if (length(marker_group_line) == 0) {
            break
        }

        marker_group_line_list <- strsplit(marker_group_line, split = c(" +", "\t"))[[1]]

        if (marker_group_line_list[1] == gene_name) {
            if (marker_group_line_list[2] == "var") {
                var <- marker_group_line_list
            } else {
                anno <- marker_group_line_list
            }
        }
    }

    lof_idx <- which(anno == "lof")
    mis_idx <- which(anno == "missense")
    syn_idx <- which(anno == "synonymous")

    lof_var <- var[lof_idx]
    mis_var <- var[mis_idx]
    syn_var <- var[syn_idx]

    out <- list(lof_var, mis_var, syn_var)
    return(out)
}


# Convert genotypes for missing / alt-first PLINK files
convert_missing <- function(plink_matrix, ref_first=FALSE) {
    if (ref_first == FALSE) {   # alt-first
        # Convert missing genotypes to reference-homozygotes
        plink_matrix[which(plink_matrix == -9)] <- 2
        # Flip genotypes
        plink_matrix <- 2 - plink_matrix
    } else {
        plink_matrix[which(plink_matrix == -9)] <- 0
    }

    return(plink_matrix)
}


# function for counting MAF of input genotype matrix (used in split_plink_matrix)
count_maf <- function(geno_mat) {
    mac_vec <- rep(0, ncol(geno_mat))
    for (i in seq_len(ncol(geno_mat))) {
        # n0 <- length(which(geno_mat[, i]) == 0)
        n1 <- length(which(geno_mat[, i] == 1))
        n2 <- length(which(geno_mat[, i] == 2))

        mac <- n1 + 2 * n2
        mac_vec[i] <- mac
    }
    maf_vec <- mac_vec / (2 * nrow(geno_mat))
    common_list <- which(maf_vec > 0.01)

    out <- list(mac_vec, common_list)
    return(out)
}


# Split genotype matrix by function annotation (lof / mis / syn)
split_plink_matrix <- function(plink_matrix, lof_var, mis_var, syn_var, mac_threshold = 10) {
    lof_mat <- plink_matrix[, lof_var]
    mis_mat <- plink_matrix[, mis_var]
    # maf_lof <- count_maf(lof_mat)[[2]]
    syn_mat <- plink_matrix[, syn_var]

    # count MAF for each functional annotation
    mac_lof <- count_maf(lof_mat)[[1]]
    common_list_lof <- count_maf(lof_mat)[[2]]

    mac_mis <- count_maf(mis_mat)[[1]]
    # maf_mis <- count_maf(mis_mat)[[2]]
    common_list_mis <- count_maf(mis_mat)[[2]]

    mac_syn <- count_maf(syn_mat)[[1]]
    # maf_syn <- count_maf(syn_mat)[[2]]
    common_list_syn <- count_maf(syn_mat)[[2]]

    # collapsing columns with MAC < threshold
    collapse_lof <- which(mac_lof < mac_threshold)
    collapse_mis <- which(mac_mis < mac_threshold)
    collapse_syn <- which(mac_syn < mac_threshold)

    rowsum_lof <- rowSums(lof_mat[, collapse_lof])
    rowsum_mis <- rowSums(mis_mat[, collapse_mis])
    rowsum_syn <- rowSums(syn_mat[, collapse_syn])

    rowsum_lof[which(rowsum_lof > 1)] <- 1
    rowsum_mis[which(rowsum_mis > 1)] <- 1
    rowsum_syn[which(rowsum_syn > 1)] <- 1

    # genotype matrix after collapsing
    lof_mat_collapsed <- cbind(lof_mat[, -c(collapse_lof, common_list_lof)], rowsum_lof)
    mis_mat_collapsed <- cbind(mis_mat[, -c(collapse_mis, common_list_mis)], rowsum_mis)
    syn_mat_collapsed <- cbind(syn_mat[, -c(collapse_syn, common_list_syn)], rowsum_syn)

    out <- list(lof_mat_collapsed, mis_mat_collapsed, syn_mat_collapsed)
    return(out)
}




# Read phenotype file
read_pheno <- function(pheno_file, pheno_code, iid_col = "f.eid") {
    pheno <- data.table::fread(pheno_file, quote = "")
    out <- subset(pheno, select = c(iid_col, pheno_code))
    return(out)
}


calc_log_lik <- function(delta, S, UtY, Y_UUtY) {
    k <- length(S)
    n <- nrow(Y_UUtY)

    log_lik1 <- 0
    for (i in 1:k) {
        log_lik1 <- log_lik1 + (UtY[i, ])^2 / (S[i] + delta)
    }

    log_lik2 <- 1 / delta * sum((Y_UUtY)^2)
    
    out <- -0.5 * (n * log(2 * pi) + sum(log(S + delta)) + (n - k) * log(delta)
                   + n + n * log(1 / n * (log_lik1 + log_lik2)))

    return(as.numeric(out))
}

calc_post_beta <- function(K, G, delta, S, UtY, U) {
    K_sparse <- as(K, "dgCMatrix")

    out <- K_sparse %*% t(G) %*% U %*% diag(1 / (S + delta)) %*% (UtY)
    return(out)
}

# Run FaST-LMM to obtain posterior beta
fast_lmm <- function(G, Y) {
    Y <- as.matrix(Y)
    K <- diag(1, nrow = ncol(G))
    L <- chol(K)
    W <- G %*% L
    W_sparse <- as(W, "dgCMatrix")
    svd_mat <- sparsesvd::sparsesvd(W_sparse)

    U <- svd_mat$u
    S <- (svd_mat$d)^2

    UtY <- t(U) %*% Y

    # Y_UUtY = Y - UUtY
    Y_UUtY <- Y - U %*% UtY

    opt <- optim(par = 1, fn = calc_log_lik, S = S, UtY = UtY, Y_UUtY = Y_UUtY,
                 method = c("Brent"), lower = 0, upper = 10000, control = list(fnscale = -1))
    opt_delta <- opt$par

    post_beta <- calc_post_beta(K, G, opt_delta, S, UtY, U)
    
    return (post_beta)
}

calc_gene_effect_size <- function(G, lof_ncol, post_beta, Y) {
    beta_lof <- post_beta[1:lof_ncol]
    beta_lof <- abs(beta_lof)
    lof_prs <- G[, 1:lof_ncol, drop = F] %*% beta_lof
    lof_prs_norm <- (lof_prs - mean(lof_prs)) / sd(lof_prs)
    m1 <- lm(Y ~ lof_prs_norm[, 1])

    return(m1$coefficients[2])
}
