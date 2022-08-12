# Description: preprocess data for running LDSC correlations across kidney traits.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(tidyverse)
library(data.table)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/hematuria_project"

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
    RSID,
    CHR,
    BP = pos,
    A1 = Allele1,
    A2 = Allele2,
    BETA = Effect,
    SE = StdErr,
    P = P.value,
    N = n)

write.table(eGFR_crea, here(project_dir, "MR_analysis", "ldsc_corr", "stanzick2021_EUR_eGFRcrea.tsv"), sep = "\t", quote = F, row.names = F)

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
    RSID,
    CHR,
    BP = pos,
    A1 = Allele1,
    A2 = Allele2,
    BETA = Effect,
    SE = StdErr,
    P = P.value,
    N = n)

write.table(eGFR_cys, here(project_dir, "MR_analysis", "ldsc_corr", "stanzick2021_EUR_eGFRcys.tsv"), sep = "\t", quote = F, row.names = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                       CKDGEN EUR ACR no DM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "MR_analysis/CKDGEN_sumstats/ACR_CKDGEN_sumstats_nonDM.csv")

ACR_noDM <- fread(file_path) %>%
  dplyr::mutate_at(c("allele1","allele2"), stringr::str_replace, "a", "A") %>%
  dplyr::mutate_at(c("allele1","allele2"), stringr::str_replace, "g", "G") %>%
  dplyr::mutate_at(c("allele1","allele2"), stringr::str_replace, "c", "C") %>%
  dplyr::mutate_at(c("allele1","allele2"), stringr::str_replace, "t", "T") %>%
  dplyr::select(
    RSID = rsID,
    A1 = allele1, # effect allele
    A2 = allele2,
    BETA = beta,
    SE = se,
    P = pval,
    N)

write.table(ACR_noDM, here(project_dir, "MR_analysis", "ldsc_corr", "CKDGEN_EUR_ACR_noDM.tsv"), sep = "\t", quote = F, row.names = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                       CKDGEN EUR ACR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "MR_analysis/CKDGEN_sumstats/ACR_CKDGEN_sumstats_all.csv")

ACR <- fread(file_path) %>%
  dplyr::mutate_at(c("allele1","allele2"), stringr::str_replace, "a", "A") %>%
  dplyr::mutate_at(c("allele1","allele2"), stringr::str_replace, "g", "G") %>%
  dplyr::mutate_at(c("allele1","allele2"), stringr::str_replace, "c", "C") %>%
  dplyr::mutate_at(c("allele1","allele2"), stringr::str_replace, "t", "T") %>%
  dplyr::select(
    RSID = rsID,
    A1 = allele1, # effect allele
    A2 = allele2,
    BETA = beta,
    SE = se,
    P = pval,
    N)

write.table(ACR, here(project_dir, "MR_analysis", "ldsc_corr", "CKDGEN_EUR_ACR.tsv"), sep = "\t", quote = F, row.names = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                       CKDGEN EUR eGFRcys
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "MR_analysis/CKDGEN_sumstats/eGFRcys_CKDGEN_sumstats.csv")

eGFRcys_ckdgen <- fread(file_path) %>%
  dplyr::mutate_at(c("allele1","allele2"), stringr::str_replace, "a", "A") %>%
  dplyr::mutate_at(c("allele1","allele2"), stringr::str_replace, "g", "G") %>%
  dplyr::mutate_at(c("allele1","allele2"), stringr::str_replace, "c", "C") %>%
  dplyr::mutate_at(c("allele1","allele2"), stringr::str_replace, "t", "T") %>%
  dplyr::select(
    RSID = rsID,
    A1 = allele1, # effect allele
    A2 = allele2,
    BETA = beta,
    SE = se,
    P = pval,
    N)

write.table(eGFRcys_ckdgen, here(project_dir, "MR_analysis", "ldsc_corr", "CKDGEN_EUR_eGFRcys.tsv"), sep = "\t", quote = F, row.names = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                       UKBB ACR      
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "MR_analysis/ACR_sumstats/ACR_GWAS_UKBB.txt")

ACR_UKBB <- fread(file_path) %>%
  mutate(CHR = as.factor(chrom), N = 395906, P = 10^(-neg_log_pvalue)) %>%
  dplyr::select(
    CHR,
    BP = pos,
    RSID = rsid,
    A1 = alt, # effect allele
    A2 = ref,
    BETA = beta,
    SE = stderr_beta,
    P,
    N)

write.table(ACR_UKBB, here(project_dir, "MR_analysis", "ldsc_corr", "UKBB_EUR_ACR.tsv"), sep = "\t", quote = F, row.names = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                       UKBB composite
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

chr_files <- list.files(here(project_dir, "MR_analysis","composite"),pattern = "^ratio.casecontrol.ukb_", full.names = T, all.files = T)

composite_list <- list()

for i in 1:length(chr_files) {

composite_list[[i]] <- fread(file_path) %>%
  mutate(CHR = as.factor(CHR)) %>%
  dplyr::select(
    CHR,
    BP = POS,
    RSID = rsid,
    A1 = Allele2, # effect allele
    A2 = Allele1,
    BETA,
    SE,
    P = p.value,
    N)

}

composite <- rbindlist(composite_list)

write.table(composite, here(project_dir, "MR_analysis", "ldsc_corr", "UKBB_EUR_composite.tsv"), sep = "\t", quote = F, row.names = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                       UKBB hematuria
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "summary_statistics", "hematuria_sumstats.tsv")

hematuria <- fread(file_path) %>%
  mutate(CHR = as.factor(chrom), N = 16235 + 378356) %>%
  dplyr::select(
    CHR,
    BP = pos,
    RSID = rsids,
    A1 = alt, # effect allele
    A2 = ref,
    BETA = beta,
    SE = sebeta,
    P = pval,
    N)

write.table(hematuria, here(project_dir, "MR_analysis", "ldsc_corr", "UKBB_EUR_hematuria.tsv"), sep = "\t", quote = F, row.names = F)
