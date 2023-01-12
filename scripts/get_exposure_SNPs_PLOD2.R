# Description - obtain exposures for the P3H2
# Selected instruments based on eQTLs from GTEx (selected based on most
# significant tissue in GTEx web browser and using all independent predictors in that tissue for
# that gene in the mashr models for TWAS in MetaXcan)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load libraries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(tidyverse)
library(stringr)
library(data.table)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir <- here("/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project")
snps_main <- c("chr3_146150702_G_T_b38","chr3_146155352_A_T_b38","chr3_146157262_G_A_b38")
snps_proxy_palindromic <- c("chr3_146140753_A_G_b38")
snps_proxy_ckdgen <- c("chr3_146073016_G_T_b38","chr3_146114877_G_A_b38")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Skeletal Muscle eQTLs on chromosome 3 only on EUR inds.
eqtls <- fread(here(project_dir,"MR_analysis","GTEx_data","Muscle_Skeletal.v8.EUR.allpairs.chr3.csv"))

# SNP list with rsids and variant_IDs
rsids <- read.table(here(project_dir,"MR_analysis","SNPs_list_PLOD2.txt"), sep = "\t", header = T)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract SNPs of interest and get rsid
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

list_all_snps <- c(snps_main, snps_proxy_palindromic, snps_proxy_ckdgen)

eqtls <- eqtls %>%
  filter(., phenotype_id == "ENSG00000152952.11") %>%
  filter(., variant_id %in% list_all_snps) %>%
  left_join(., rsids, by = "variant_id")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Format tables to use in MRbase analysis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# The effect allele is the second allele in the variant ID
# Need to find a way to systematically include rsid instead of variant_ID

eqtls <- eqtls %>%
  mutate(Phenotype = "PLOD2_eQTL_Skeletal_Muscle",
         alleles = stringr::str_extract(variant_id, "[:upper:][:punct:][:upper:]")) %>%
  separate(alleles, c("other_allele","effect_allele"),sep="_") %>%
  select(
    Phenotype,
    SNP = rsid,
    effect_allele,
    other_allele,
    beta = slope,
    se = slope_se,
    pval = pval_nominal
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write outcome files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.table(eqtls, here(project_dir,"MR_analysis","MR_exposures","Skeletal_Muscle_EUR_PLOD2_exposure.tsv"),
	sep = "\t", row.names = F, quote = F)



