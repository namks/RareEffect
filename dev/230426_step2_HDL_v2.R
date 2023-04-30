# Argument = chromosome
args <- commandArgs(trailingOnly = TRUE)
chr <- args[1]
library(data.table, quietly = TRUE)
library(seqminer, quietly = TRUE)
source("/media/leelabsg-storage0/kisung/RVPRS/dev/src/RVPRS_function.R")

# Read SAIGE-GENE+ results (region-level summary statistics)
setwd("/media/leelabsg-storage0/kisung/dnanexus/WES_470k/output/conti")
gene_result <- paste0("SAIGE_GENE_WES_470k_f.30760.0.0_chr", chr, "_WB_withQC.txt")
d <- fread(gene_result)

# Load step 1 results
step1_rda_name <- "/media/leelabsg-storage0/kisung/RVPRS/dev/output/step1/step1_HDL.rda"
load(step1_rda_name)
sigma_sq <- var(modglmm$residuals)
n_samples <- length(modglmm$residuals)

# Filter genes/regions with p-value cutoff
filter_groupname <- "Cauchy"
pval_cutoff <- 0.05
mac_threshold <- 10

genes <- d[which((d$Group == filter_groupname) & (d$Pvalue < pval_cutoff)),]$Region

# Dataframe for outputs
single_var_output <- NULL
single_outname <- paste0("/media/leelabsg-storage0/kisung/RVPRS/dev/output/step2/step2_HDL_chr", chr, "_single.txt")
gene_output <- NULL
gene_outname <- paste0("/media/leelabsg-storage0/kisung/RVPRS/dev/output/step2/step2_HDL_chr", chr, "_gene.txt")

# Iterate by gene/region
for (g in genes) {
    time_outname <- paste0("/media/leelabsg-storage0/kisung/RVPRS/dev/output/time/step2_HDL_chr", chr, "_", g, ".time")
    print(paste0("Gene ", g, " started"))
    # if (file.exists(single_outname)) {
    #     next
    # }
    time_elapsed <- system.time({
        # Read groupfile to read range and functional annotations
        groupfile_name <- "/media/leelabsg-storage0/kisung/dnanexus/group_files/UKBexome_all_chr.txt"

        var_by_func_anno <- read_groupfile(groupfile_name, g)
        if (length(var_by_func_anno[[1]]) == 0)
            next
        pos_range <- get_range(var_by_func_anno)

        plink_prefix <- paste0("/media/leelabsg-storage0/DATA/UKBB/WES/UKBB_WES_200k/ukb23155_c", chr, "_b0_v1")
        bim <- fread(paste0(plink_prefix, ".bim"))
        fam <- fread(paste0(plink_prefix, ".fam"))
        start_idx <- min(which(bim$V4 == pos_range[[1]]))
        end_idx <- max(which(bim$V4 == pos_range[[2]]))

        # Need to revise load_plink function in RVPRS_function.R
        plink_matrix <- seqminer::readPlinkToMatrixByIndex(plink_prefix, sampleIndex = 1:nrow(fam), markerIndex = start_idx:end_idx)
        colnames(plink_matrix) <- bim$V2[start_idx:end_idx]
        plink_matrix_converted <- convert_missing(plink_matrix, ref_first = FALSE)
        mat_by_func_anno <- split_plink_matrix(plink_matrix_converted,
                            var_by_func_anno[[1]],
                            var_by_func_anno[[2]],
                            var_by_func_anno[[3]],
                            mac_threshold = mac_threshold
                        )

        # Make genotype matrix
        G <- cbind(mat_by_func_anno[[1]], mat_by_func_anno[[2]], mat_by_func_anno[[3]])
        lof_ncol <- ncol(mat_by_func_anno[[1]])
        mis_ncol <- ncol(mat_by_func_anno[[2]])
        syn_ncol <- ncol(mat_by_func_anno[[3]])

        y_tilde <- cbind(as.numeric(modglmm$sampleID), modglmm$residuals)

        G_reordered <- G[match(y_tilde[,1], rownames(G)),]

        post_beta_lof <- NULL
        post_beta_mis <- NULL
        post_beta_syn <- NULL

        # Refactorize this code later -> make function in RVPRS_function.R
        if (lof_ncol > 0) {
            fast_lmm_lof <- fast_lmm(G = G_reordered[,c(1:lof_ncol), drop = F], Y = y_tilde[,2])
            post_beta_lof <- fast_lmm_lof[[1]]
            tr_GtG <- fast_lmm_lof[[2]]
            delta <- fast_lmm_lof[[3]]
            h2_lof <- (sigma_sq / delta) * tr_GtG / ((sigma_sq / delta) * tr_GtG + n_samples * sigma_sq)
        }
        if (mis_ncol > 0) {
            fast_lmm_mis <- fast_lmm(G = G_reordered[,c((lof_ncol + 1):(lof_ncol + mis_ncol)), drop = F], Y = y_tilde[,2])
            post_beta_mis <- fast_lmm_mis[[1]]
            tr_GtG <- fast_lmm_mis[[2]]
            delta <- fast_lmm_mis[[3]]
            h2_mis <- (sigma_sq / delta) * tr_GtG / ((sigma_sq / delta) * tr_GtG + n_samples * sigma_sq)
        }
        if (syn_ncol > 0) {
            fast_lmm_syn <- fast_lmm(G = G_reordered[,c((lof_ncol + mis_ncol + 1):(lof_ncol + mis_ncol + syn_ncol)), drop = F], Y = y_tilde[,2])
            post_beta_syn <- fast_lmm_syn[[1]]
            tr_GtG <- fast_lmm_syn[[2]]
            delta <- fast_lmm_syn[[3]]
            h2_syn <- (sigma_sq / delta) * tr_GtG / ((sigma_sq / delta) * tr_GtG + n_samples * sigma_sq)
        }
        
        # output related to single-variant effect size
        post_beta <- rbind(post_beta_lof, post_beta_mis, post_beta_syn)
        func_anno <- c(rep("lof", lof_ncol), rep("mis", mis_ncol), rep("syn", syn_ncol))
        single_effect <- cbind(chr, g, func_anno, colnames(G), post_beta[,1])
        colnames(single_effect) <- c("CHR", "REGION", "FUNC_ANNO", "SNP", "BETA")

        # output related to gene-level heritability
        h2_gene <- c(h2_lof, h2_mis, h2_syn, h2_lof + h2_mis + h2_syn)
        gene_effect <- cbind(chr, g, c("lof", "mis", "syn", "sum"), h2_gene)
        colnames(gene_effect) <- c("CHR", "REGION", "FUNC_ANNO", "h2")
    })

    write.table(t(time_elapsed), time_outname, row.names=F, col.names=T, quote=F)
    # write.table(single_effect, single_outname, row.names=F, col.names=F, quote=F)
    single_var_output <- rbind(single_var_output, single_effect)
    gene_output <- rbind(gene_output, gene_effect)
    print(paste0("Gene ", g, " completed."))
}

write.table(single_var_output, single_outname, row.names=F, col.names=T, quote=F)
write.table(gene_output, gene_outname, row.names=F, col.names=T, quote=F)
print(paste0("Chromosome ", chr, " completed."))
