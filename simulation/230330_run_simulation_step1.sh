#!/bin/bash

for n in {1..100}
do
  for i in {1..2}
  do
    for j in {1..2}
    do
      for k in {1..2}
      do
        # execute the docker command with n, i, j, and k as arguments
        docker run -v /media/leelabsg-storage0:/media/leelabsg-storage0 wzhou88/saige:1.1.3 step1_fitNULLGLMM.R \
        --plinkFile=/media/leelabsg-storage0/kisung/RVPRS/simulation/common_variants/common_variants_WB \
        --phenoFile=/media/leelabsg-storage0/kisung/RVPRS/simulation/data/pheno/merged/pheno_sim${n}_${i}_${j}_${k}.txt \
        --phenoCol=pheno_total \
        --covarColList=cov1,cov2 \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative \
        --invNormalize=TRUE \
        --outputPrefix=/media/leelabsg-storage0/kisung/RVPRS/simulation/script/result/step1/step1_sim${n}_${i}_${j}_${k} \
        --nThreads=8 \
        --LOCO=FALSE \
        --IsOverwriteVarianceRatioFile=TRUE
      done
    done
  done
done
