# Description: prepare data for colocalization analysis across kidney disease traits on genes PLOD1 and PLOD2

# 1. PLOD1 expression and eGFRcrea (meta-analysis)
# 2. PLOD1 expression and eGFRcys (meta-analysis)
# 3. PLOD2 expression and eGFRcys (CKDGen)

# Packages -------------------------------------------------------

library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(arrow)
library(colochelpR)
library(rutils)
library(BiocManager)
library(data.table)
library(biomaRt)
library(BSgenome)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)

# Arguments ----------------------------------------------------------------------

project_dir <- "/home/fridald4/projects/def-gsarah/fridald4/renal_genetics_project/MR_analysis"

dbsnp_144.GRCh37 <- SNPlocs.Hsapiens.dbSNP144.GRCh37
dbsnp_144.GRCh38 <- SNPlocs.Hsapiens.dbSNP144.GRCh38

# PLOD1 coordinates (GRCh37)
plod1_low_lim <- 11994774 - 1000000
plod1_upp_lim <- 12035594 + 1000000

# PLOD1 coordinates (GRCh38)
plod1_low_lim2 <- 11934717 - 1000000
plod1_upp_lim2 <- 11975537 + 1000000

# PLOD2 coordinates (GRCh37)
plod2_low_lim <- 145787227 - 1000000
plod2_upp_lim <- 145878971 + 1000000

# PLOD2 coordinates (GRCh38)
plod2_low_lim2 <- 146069440 - 1000000
plod2_upp_lim2 <- 146161184 + 1000000

# P3H2 coordinates (GRCh37)
p3h2_low_lim <- 189674517 - 1000000
p3h2_upp_lim <- 189840067 + 1000000

# P3H2 coordinates (GRCh38)
p3h2_low_lim2 <- 189956728 - 1000000
p3h2_upp_lim2 <- 190122278 + 1000000

# Read and format data files -------------------------------------------------------

# Note: for all GWAS/QTL traits will use RSIDs in colocalization (not chr:pos), to avoid issues with genome build versions.

# GWAS summary statistics
eGFRcrea_meta <- read.table(here(project_dir, "stanzick2021_meta_analysis_EUR_eGFRcrea_CKDGEN_UKBB.gc"), sep = "\t", header = T) %>%
  mutate(GWAS = "eGFRcrea_stanzick2021",
         maf = case_when(
           Freq1 <= 0.5 ~ Freq1,
           Freq1 > 0.5 ~ 1-Freq1)) %>%
  dplyr::select(GWAS, SNP = RSID, CHR = chr, BP = pos, beta = Effect, se = StdErr, p.value = P.value, A1 = Allele1, A2 = Allele2, maf, N = n) %>%
  colochelpR::get_varbeta(.)

head(eGFRcrea_meta)

eGFRcys_meta <- read.table(here(project_dir, "stanzick2021_meta_analysis_EUR_eGFRcys_CKDGEN_UKBB.gc"), sep = "\t", header = T)  %>%
  mutate(GWAS = "eGFRcys_stanzick2021",
         maf = case_when(
           Freq1 <= 0.5 ~ Freq1,
           Freq1 > 0.5 ~ 1-Freq1)) %>%
  dplyr::select(GWAS, SNP = RSID, CHR = chr, BP = pos, beta = Effect, se = StdErr, p.value = P.value, A1 = Allele1, A2 = Allele2, maf, N = n) %>%
  colochelpR::get_varbeta(.)

head(eGFRcys_meta)

eGFRcys_ckdgen <- read.table(here(project_dir, "CKDGEN_sumstats", "eGFRcys_CKDGEN_sumstats.csv"), sep = ",", header = T) %>%
  colochelpR::convert_rs_to_loc(df = ., SNP_column = "rsID", dbsnp_144.GRCh37) %>%
  mutate(GWAS = "eGFRcys_CKDGEN",
         maf = case_when(
           freqA1 <= 0.5 ~ freqA1,
           freqA1 > 0.5 ~ 1-freqA1)) %>%
  separate(loc, c("CHR", "BP"), sep = ":") %>%
  dplyr::select(GWAS, SNP = rsID, CHR, BP, beta, se, p.value = pval, A1 = allele1, A2 = allele2, maf, N) %>%
  colochelpR::get_varbeta(.)

head(eGFRcys_ckdgen)

hematuria_ukbb <- read.table("/home/fridald4/projects/def-gsarah/fridald4/renal_genetics_project/summary_statistics/hematuria_sumstats.tsv", sep = "\t", header = T) %>%
  mutate(GWAS = "hematuria_UKBB",
         N = 16235 + 378356,
         maf = case_when(
           af <= 0.5 ~ af,
           af > 0.5 ~ 1-af)) %>%
  dplyr::select(GWAS, SNP = rsids, CHR = chrom, BP = pos, beta, se = sebeta, p.value = pval, A1 = alt, A2 = ref, maf, N) %>%
  colochelpR::get_varbeta(.)

# PLOD1 (ENSG00000083444.16 - chr1) GTEx eQTLs (whole blood) GRCh38
# GRCh37 coordinates: 145787227..145878971
# GRCh38 coordinates: 146069440..146161184
plod1 <- read_parquet(here(project_dir, "GTEx_data", "Whole_Blood.v8.EUR.allpairs.chr1.parquet")) %>%
  # keep only PLOD1 eQTLs
  filter(., str_detect(phenotype_id, "ENSG00000083444.16")) %>%
  separate(., variant_id, c("CHR","BP","NEA","EA",NA), sep = "_") %>%
  mutate(., CHR = as.integer(str_remove_all(CHR, "chr")), eQTL_dataset = "PLOD1_expression", N = 573) %>%
  # get rsids
  colochelpR::convert_loc_to_rs(df = ., dbsnp_144.GRCh38)

  # remove duplicate rsids
  plod1 <- plod1 %>% 
  dplyr::mutate(CHR_BP = stringr::str_c(CHR, ":", BP)) %>%
  dplyr::group_by(CHR_BP) %>%
  dplyr::filter(!any(row_number() > 1))  %>%
  dplyr::ungroup() %>%
  dplyr::select(
    eQTL_dataset,
    gene = phenotype_id,
    SNP,
    CHR,
    BP,
    beta = slope,
    se = slope_se,
    p.value = pval_nominal,
    A1 = EA,
    A2 = NEA,
    maf,
    N
  ) %>%
  colochelpR::get_varbeta(.)

head(plod1)

# PLOD2 (ENSG00000152952.11 - chr3) GTEx eQTLs (skeletal muscle) GRCh38
# GRCh37 coordinates: 11994774..12035594
# GRCh38 coordinates: 11934717..11975537
plod2 <- read_parquet(here(project_dir, "GTEx_data", "Muscle_Skeletal.v8.EUR.allpairs.chr3.parquet")) %>%
  # keep only PLOD2 eQTLs
  filter(., str_detect(phenotype_id, "ENSG00000152952.11")) %>%
  separate(., variant_id, c("CHR","BP","NEA","EA",NA), sep = "_") %>%
  mutate(., CHR = as.integer(str_remove_all(CHR, "chr")), eQTL_dataset = "PLOD2_expression", N = 602) %>%
  # get rsids
  colochelpR::convert_loc_to_rs(df = ., dbsnp_144.GRCh38)

  # remove duplicate rsids
  plod2 <- plod2 %>% 
  dplyr::mutate(CHR_BP = stringr::str_c(CHR, ":", BP)) %>%
  dplyr::group_by(CHR_BP) %>%
  dplyr::filter(!any(row_number() > 1))  %>%
  dplyr::ungroup() %>%
  dplyr::select(
    eQTL_dataset,
    gene = phenotype_id,
    SNP,
    CHR,
    BP,
    beta = slope,
    se = slope_se,
    p.value = pval_nominal,
    A1 = EA,
    A2 = NEA,
    maf,
    N
  ) %>%
  colochelpR::get_varbeta(.)

head(plod2)

# P3H2 (ENSG00000090530.9 - chr3) GTEx eQTLs (Nerve Tibial) GRCh38
# GRCh37 coordinates: 189674517..189840067
# GRCh38 coordinates: 189956728..190122278

p3h2 <- read_parquet(here(project_dir, "GTEx_data", "Nerve_Tibial.v8.EUR.allpairs.chr3.parquet")) %>%
  # keep only P3H2 eQTLs
  filter(., str_detect(phenotype_id, "ENSG00000090530.9")) %>%
  separate(., variant_id, c("CHR","BP","NEA","EA",NA), sep = "_") %>%
  mutate(., CHR = as.integer(str_remove_all(CHR, "chr")), eQTL_dataset = "P3H2_expression", N = 449) %>%
  # get rsids
  colochelpR::convert_loc_to_rs(df = ., dbsnp_144.GRCh38)

# remove duplicate rsids
  p3h2 <- p3h2 %>% 
  dplyr::mutate(CHR_BP = stringr::str_c(CHR, ":", BP)) %>%
  dplyr::group_by(CHR_BP) %>%
  dplyr::filter(!any(row_number() > 1))  %>%
  dplyr::ungroup() %>%
  dplyr::select(
    eQTL_dataset,
    gene = phenotype_id,
    SNP,
    CHR,
    BP,
    beta = slope,
    se = slope_se,
    p.value = pval_nominal,
    A1 = EA,
    A2 = NEA,
    maf,
    N
  ) %>%
  colochelpR::get_varbeta(.)

head(p3h2)

# Get 1Mb regions flanking the genes -------------------------------------------------------

eGFRcrea_meta_PLOD1 <- eGFRcrea_meta %>%
  filter(., CHR == 1) %>%
  filter(., between(BP, plod1_low_lim, plod1_upp_lim))

cat("The number of SNPs used in eGFRcrea_meta for PLOD1 are: ", nrow(eGFRcrea_meta_PLOD1), ".\n")

eGFRcys_meta_PLOD1 <- eGFRcys_meta %>%
  filter(., CHR	== 1) %>%
  filter(., between(BP, plod1_low_lim, plod1_upp_lim))

cat("The number of SNPs used in eGFRcys_meta for PLOD1 are: ", nrow(eGFRcys_meta_PLOD1), ".\n")

#plod1_qtls <- plod1 %>%
#  filter(., between(BP, plod1_low_lim2, plod1_upp_lim2))

plod1_qtls <- plod1 %>%
  filter(., SNP %in% unique(c(eGFRcrea_meta_PLOD1$SNP, eGFRcys_meta_PLOD1$SNP)))

cat("The number of QTLs for PLOD1 are: ", nrow(plod1_qtls), ".\n")

eGFRcys_ckdgen_PLOD2 <- eGFRcys_ckdgen %>%
  filter(., CHR	== 3) %>%
  filter(., between(BP, plod2_low_lim, plod2_upp_lim))

cat("The number of SNPs used in eGFRcys_CKDGEN for PLOD2 are: ", nrow(eGFRcys_ckdgen_PLOD2), ".\n")

# plod2_qtls <- plod2 %>%
#  filter(., between(BP, plod2_low_lim2, plod2_upp_lim2))

plod2_qtls <- plod2 %>%
  filter(., SNP %in% eGFRcys_ckdgen_PLOD2$SNP)

cat("The number of QTLs for PLOD2 are: ", nrow(plod2_qtls), ".\n")

hematuria_ukbb_P3H2 <- hematuria_ukbb %>%
  filter(., CHR == 3) %>%
  filter(., between(BP, p3h2_low_lim2, p3h2_upp_lim2))

cat("The number of SNPs used in hematuria_ukbb for P3H2 are: ", nrow(hematuria_ukbb_P3H2), ".\n")

p3h2_qtls <- p3h2 %>%
  filter(., SNP %in% hematuria_ukbb_P3H2$SNP)

# Save files for coloc -------------------------------------------------------

write.table(eGFRcrea_meta_PLOD1, here(project_dir, "colocalization", "data_for_coloc", "eGFRcrea_meta_PLOD1.tsv"),
            sep = "\t", row.names = F, quote = F)
write.table(eGFRcys_meta_PLOD1, here(project_dir, "colocalization", "data_for_coloc", "eGFRcys_meta_PLOD1.tsv"),
            sep = "\t", row.names = F, quote = F)
write.table(plod1_qtls, here(project_dir, "colocalization", "data_for_coloc", "PLOD1_qtls.tsv"),
            sep = "\t", row.names = F, quote = F)
write.table(eGFRcys_ckdgen_PLOD2, here(project_dir, "colocalization", "data_for_coloc", "eGFRcys_meta_PLOD2.tsv"),
            sep = "\t", row.names = F, quote = F)
write.table(plod2_qtls, here(project_dir, "colocalization", "data_for_coloc", "PLOD2_qtls.tsv"),
            sep = "\t", row.names = F, quote = F)
write.table(hematuria_ukbb_P3H2, here(project_dir, "colocalization", "data_for_coloc", "hematuria_ukbb_P3H2.tsv"),
            sep = "\t", row.names = F, quote = F)
write.table(p3h2_qtls, here(project_dir, "colocalization", "data_for_coloc", "P3H2_qtls.tsv"),
            sep = "\t", row.names = F, quote = F)
