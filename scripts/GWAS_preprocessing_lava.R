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
# library(devtools)
# library(BiocManager)
library(dplyr)
library(tidyr)
library(data.table)
# library(biomaRt)
# library(BSgenome)
# library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

# devtools::install_github("RHReynolds/colochelpR")
# devtools::install_github("RHReynolds/rutils")
# library(colochelpR)
# library(rutils)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/renal_genetics_project"
# dbsnp_144 <- SNPlocs.Hsapiens.dbSNP144.GRCh37
# path_hg38_hg37 <- here(project_dir,"reference_data","hg38ToHg19.over.chain")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#					eGFR creatinine
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "MR_analysis/stanzick2021_meta_analysis_EUR_eGFRcrea_CKDGEN_UKBB.gc")

eGFR_crea <- fread(file_path) %>%
  dplyr::mutate(CHR = as.factor(chr))  %>%
  dplyr::mutate_at(c("Allele1","Allele2"), stringr::str_replace, "a", "A") %>%
  dplyr::mutate_at(c("Allele1","Allele2"), stringr::str_replace, "g", "G") %>%
  dplyr::mutate_at(c("Allele1","Allele2"), stringr::str_replace, "c", "C") %>%
  dplyr::mutate_at(c("Allele1","Allele2"), stringr::str_replace, "t", "T") %>%
  dplyr::select(
    SNP = RSID,
    CHR,
    BP = pos,
    A1 = Allele1,
    A2 = Allele2,
    BETA = Effect,
    SE = StdErr,
    P = P.value,
    N = n)

write.table(eGFR_crea, here(project_dir, "MR_analysis", "lava", "input_data", "stanzick2021_EUR_eGFRcrea.lava.gz"), sep = "\t", quote = F, row.names = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                       eGFR cystatine
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "MR_analysis/stanzick2021_meta_analysis_EUR_eGFRcys_CKDGEN_UKBB.gc")

eGFR_cys <- fread(file_path) %>%
  dplyr::mutate(CHR = as.factor(chr))  %>%
  dplyr::mutate_at(c("Allele1","Allele2"), stringr::str_replace, "a", "A") %>%
  dplyr::mutate_at(c("Allele1","Allele2"), stringr::str_replace, "g", "G") %>%
  dplyr::mutate_at(c("Allele1","Allele2"), stringr::str_replace, "c", "C") %>%
  dplyr::mutate_at(c("Allele1","Allele2"), stringr::str_replace, "t", "T") %>%
  dplyr::select(
    SNP = RSID,
    CHR,
    BP = pos,
    A1 = Allele1,
    A2 = Allele2,
    BETA = Effect,
    SE = StdErr,
    P = P.value,
    N = n)

write.table(eGFR_cys, here(project_dir, "MR_analysis", "lava", "input_data", "stanzick2021_EUR_eGFRcys.lava.gz"), sep = "\t", quote = F, row.names = F)
