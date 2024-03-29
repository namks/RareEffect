# RareEffect

This is a repository for the test code of RareEffect

## Simulation study

### Phenotype Generation Scheme

1. Phenotype was generated using real genotype data from UK Biobank (UKB).
2. Used randomly selected 30,000 SNPs after LD pruning from UKB called genotype data.
3. For rare variants, used 10 genes: "ANGPTL3", "ANGPTL8", "APOA1", "CETP", "DOCK6", "FAM86C1", "LPL", "PLA2G12A", "PLIN1", "SMARCA4"
4. Generated phenotypes under generative model from SAIGE-GENE+ using 8 (2 * 2 * 2) scenarios
5. Repeated phenotype generation process for 100 times

#### Generation process

1. Run `230329_generate_pheno_common.R` and `230329_generate_pheno_rare.R`
2. Run `230329_merge_generated_pheno.R` to merge generated phenotypes by simulation scenario

### Run Analysis

This method consists of three steps.

Step 1. Fitting the null model

1. Run `230330_run_simulation_step1.sh`

Step 2. Estimate the variant-level effect size and gene-level effect size for rare variants

1. Run `230330_run_simulation_step2.R`

#### Run Analysis for Baseline Methods

Baseline methods contains

* Linear regression (linear regression with multiple variants)
* Marginal effect (linear regression with a single variant)
* Ridge regression
