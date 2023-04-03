setwd("~/Downloads/corr/")
library(tidyr)
library(dplyr)
library(ggplot2)

corr_all <- NULL
func_anno <- c("all", "lof", "mis", "syn")

for (i in 1:10) {
    for (j in 1:4) {
        for (ii in 1:2) {
            for (jj in 1:2) {
                for (kk in 1:2) {
                    fname <- paste0("corr_", func_anno[j], "_gene_", i, "_", ii, "_", jj, "_", kk, ".txt")
                    d <- read.table(fname)
                    temp <- cbind(func_anno[j], i, ii, jj, kk, d$V1)
                    corr_all <- rbind(corr_all, temp)
                }
            }
        }
    }
}

corr_all <- as.data.frame(corr_all)
corr_all$corr <- as.numeric(corr_all$corr)
colnames(corr_all) <- c("func_anno", "gene", "ii", "jj", "kk", "corr")
head(corr_all)

corr_all2 <- corr_all[which(corr_all$gene == 9),]
head(corr_all2)
boxplot(corr ~ func_anno, data=corr_all2)

mean(corr_all2[which(corr_all2$func_anno == "lof"),]$corr, na.rm=T)
mean(corr_all2[which(corr_all2$func_anno == "mis"),]$corr, na.rm=T)
mean(corr_all2[which(corr_all2$func_anno == "syn"),]$corr, na.rm=T)
mean(corr_all2$corr, na.rm=T)

corr_all3 <- corr_all[which((corr_all$ii == 1) & (corr_all$jj == 1) & (corr_all$kk == 1)),]
boxplot(corr ~ func_anno, data=corr_all3)

mean(corr_all3[which(corr_all3$func_anno == "lof"),]$corr, na.rm=T)
mean(corr_all3[which(corr_all3$func_anno == "mis"),]$corr, na.rm=T)
mean(corr_all3[which(corr_all3$func_anno == "syn"),]$corr, na.rm=T)
mean(corr_all3$corr, na.rm=T)


gene_score <- NULL
n_sim <- 100

for (n in 1:n_sim) {
    for (i in 1:10) {
        for (ii in 1:2) {
            for (jj in 1:2) {
                for (kk in 1:2) {
                    fname <- paste0("./gene/sim", n, "_gene", i, "_", ii, "_", jj, "_", kk, ".txt")
                    d <- read.table(fname)
                    temp <- cbind(n, i, ii, jj, kk, d$V1)
                    gene_score <- rbind(gene_score, temp)
                }
            }
        }
    }    
}
colnames(gene_score) <- c("sim", "gene", "ii", "jj", "kk", "gene_effect")
head(gene_score)
gene_score <- as.data.frame(gene_score)
gene_score2 <- gene_score[which((gene_score$ii == 2) & (gene_score$jj == 2) & (gene_score$kk == 2)),]
mean(gene_score2$gene_effect)

boxplot(gene_effect ~ ii + jj, xlab="scenario" , data=gene_score)
