library(data.table)
library(dplyr)

setwd("/media/leelabsg-storage0/kisung/RVPRS/simulation/data/pheno")

n_sim <- 100

for (n in 1:n_sim) {
    pheno_common_name <- paste0("pheno_sim", n, ".txt")
    pheno_common <- fread(pheno_common_name)
    for (ii in 1:2) {
        for (jj in 1:2) {
            for (kk in 1:2) {
                for (j in 1:10) {
                    pheno_rare_name <- paste0("pheno_sim", n, "_gene", j, "_", ii, "_", jj, "_", kk, ".txt")
                    pheno_rare <- fread(pheno_rare_name)
                    
                    if (j == 1) {
                        pheno_rare_total <- pheno_rare
                    } else {
                        pheno_rare_total$pheno_rv <- pheno_rare_total$pheno_rv + pheno_rare$pheno_rv
                    }
                }
                pheno <- left_join(pheno_common, pheno_rare_total, by="IID")
                pheno$pheno_total <- pheno$pheno_common + pheno$cov1 + pheno$cov2 + pheno$noise + pheno$pheno_rv

                outname <- paste0("./merged/pheno_sim", n, "_", ii, "_", jj, "_", kk, ".txt")
                write.table(pheno, outname, row.names=F, col.names=T, quote=F)
            }
        }
    }
    print(paste0("Simulation ", n, " merged."))
}
