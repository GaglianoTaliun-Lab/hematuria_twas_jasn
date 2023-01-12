# Description: Retrieve summary statistics from CKDGen (ACR and eGFR) and hematuria-phecode-593
# for the SNPs of interest to perform meta-analyses with the same phenotypes from UKBB

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load libraries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(tidyverse)
library(stringr)
library(data.table)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir <- here("/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project")

# snps that are present in both datasets (using the proxy for the SNPs not present in CKDGEN and for the palindromic SNP)
snps_P3H2 <- c("rs11915773","rs710590","rs6783292","rs838268")
snps_PLOD1 <- c("rs873458","rs2273291")
snps_PLOD2 <- c("rs1449444","rs1170389","rs4681295")

snps_chr1 <- c(snps_PLOD1)
snps_chr3 <- c(snps_PLOD2, snps_P3H2)

snps_of_interest <- c(snps_P3H2, snps_PLOD1, snps_PLOD2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ACR UKBB GWAS summary statistics:
ACR_UKBB_chr1 <- fread(here(project_dir,"MR_analysis","ACR_sumstats","ratio.ukb_chr1_v3.SAIGE.txt"))
ACR_UKBB_chr3 <- fread(here(project_dir,"MR_analysis","ACR_sumstats","ratio.ukb_chr3_v3.SAIGE.txt"))

# CKDGEN files (read them in a list)
CKDGEN_list <- setNames(
  object = 
    list.files(path = here(project_dir,"MR_analysis","CKDGEN_sumstats"), pattern = "*.csv", full.names = T) %>%
    lapply(., function(x)fread(x)),
  nm =
    list.files(path = here(project_dir,"MR_analysis","CKDGEN_sumstats"), pattern = "*.csv", full.names = F) %>%
    stringr::str_remove(., "_sumstats") %>%
    stringr::str_remove(., ".csv")
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# For METAL, I need rsid, chr, pos, beta, se, pvalue, ref allele, non ref allele

# ACR UKBB chr1
ACR_UKBB_chr1 <- ACR_UKBB_chr1 %>%
  filter(rsid %in% snps_chr1) %>%
  select(
    SNP = rsid,
    chr = CHR,
    pos = POS,
    other_allele = Allele1,
    effect_allele = Allele2,
    pval = p.value,
    beta = BETA,
    se = SE,
    EAF = AF_Allele2
)

# ACR UKBB chr3
ACR_UKBB_chr3 <- ACR_UKBB_chr3 %>%
  filter(rsid %in% snps_chr3) %>%
  select(
    SNP = rsid,
    chr = CHR,
    pos = POS,
    other_allele = Allele1,
    effect_allele = Allele2,
    pval = p.value,
    beta = BETA,
    se = SE,
    EAF = AF_Allele2     
)

# merge both chrs data
ACR_UKBB <- rbind(ACR_UKBB_chr1, ACR_UKBB_chr3)

# ACR CKDGEN
# rsID,allele1,allele2,freqA1,beta,se,pval,N
CKDGEN_list <- CKDGEN_list %>%
  lapply(., function(x){
    filter(x, rsID %in% snps_of_interest) %>%
    select(
      SNP = rsID,
      other_allele = allele2,
      effect_allele = allele1,
      pval = pval,
      beta = beta,
      se = se,
      EAF = freqA1
    )
  })

cat("Length of CKDGEN list: ",length(CKDGEN_list),". \n")
cat(names(CKDGEN_list))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.table(ACR_UKBB, 
            here(project_dir,"MR_analysis","meta_analysis_CKDgen_UKBB","input_ACR_UKBB_P3H2_PLOD1_PLOD2.txt"), 
            sep = "\t", row.names = F, quote = F)

for (i in 1:length(CKDGEN_list)) {
  write.table(CKDGEN_list[[i]], 
              here(project_dir,"MR_analysis","meta_analysis_CKDgen_UKBB",stringr::str_c("input_",names(CKDGEN_list)[[i]],"_P3H2_PLOD1_PLOD2.txt")), 
              sep = "\t", row.names = F, quote = F)
}
  
