# process GWAS summary statistics for LAVA - to run as a normal job on compute canada

# Summary statistics must have the following columns:
# SNP / ID / SNPID_UKB/ SNPID / MarkerName / RSID / RSID_UKB: SNP IDs
# A1 / ALT: effect allele
# A2 / REF: reference allele
# N / NMISS / OBS_CT / N_analyzed: number of samples
# Z / T / STAT / Zscore: if provided, no p-values or coefficients are needed; otherwise, please provide both:
# B / BETA / OR / logOdds: effect size coefficients
# P: p-values

# install.packages("devtools")
# install.packages("BiocManager")

# BiocManager::install("biomaRt")
# BiocManager::install("BSgenome")

library(here)
library(devtools)
library(BiocManager)
library(dplyr)
library(tidyr)
library(data.table)
library(biomaRt)
library(BSgenome)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

# devtools::install_github("RHReynolds/colochelpR")
# devtools::install_github("RHReynolds/rutils")
library(colochelpR)
library(rutils)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project"
dbsnp_144 <- SNPlocs.Hsapiens.dbSNP144.GRCh37
# path_hg38_hg37 <- here(project_dir,"reference_data","hg38ToHg19.over.chain")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#					ACR UKBB
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "MR_analysis/ldsc_corr/acr.sumstats.gz")

acr <- fread(file_path) %>%
  dplyr::select(
    SNP,
    A1,
    A2,
    Z,
    N) %>%
  colochelpR::convert_rs_to_loc(df = ., SNP_column = "SNP", dbSNP = dbsnp_144) %>%
  separate(loc, c("CHR", "BP"), sep = ":")

write.table(acr, here(project_dir, "MR_analysis", "lava", "input_data", "acr.lava.gz"), sep = "\t", quote = F, row.names = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                       eGFR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "MR_analysis/ldsc_corr/regenie-egfr.sumstats.gz")

eGFR <- fread(file_path) %>%
  dplyr::select(
    SNP,
    A1,
    A2,
    Z,
    N) %>%
  colochelpR::convert_rs_to_loc(df = ., SNP_column = "SNP", dbSNP = dbsnp_144) %>%
  separate(loc, c("CHR", "BP"), sep = ":")

write.table(eGFR, here(project_dir, "MR_analysis", "lava", "input_data", "regenie-egfr.lava.gz"), sep = "\t", quote = F, row.names = F)
