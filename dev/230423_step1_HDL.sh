docker run -dv /media/leelabsg-storage0:/media/leelabsg-storage0 wzhou88/saige:1.1.3 step1_fitNULLGLMM.R \
        --plinkFile=/media/leelabsg-storage0/GWAS_tutorial/data/UKB_step1 \
        --phenoFile=/media/leelabsg-storage0/kisung/RVPRS/dev/data/pheno_HDL.txt \
        --phenoCol=HDL \
        --covarColList=Age,Sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative \
        --invNormalize=TRUE \
        --outputPrefix=/media/leelabsg-storage0/kisung/RVPRS/dev/output/step1/step1_HDL \
        --nThreads=8 \
        --LOCO=FALSE \
        --IsOverwriteVarianceRatioFile=TRUE
