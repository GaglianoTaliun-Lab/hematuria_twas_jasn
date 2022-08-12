# Description - obtain exposures for the PLOD1
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

project_dir <- here("/home/fridald4/projects/def-gsarah/fridald4/renal_genetics_project")
snps_main <- c("chr1_11986032_G_A_b38","chr1_11949405_G_A_b38")
snps_proxy_ckdgen <- c("chr1_11965296_T_C_b38")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Nerve Tibial eQTLs on chromosome 3 only on EUR inds.
eqtls <- fread(here(project_dir,"MR_analysis","GTEx_data","Whole_Blood.v8.EUR.allpairs.chr1.csv"))

# SNP list with rsids and variant_IDs
rsids <- read.table(here(project_dir,"MR_analysis","SNPs_list_PLOD1.txt"), sep = "\t", header = T)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract SNPs of interest and get rsid
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

list_all_snps <- c(snps_main, snps_proxy_ckdgen)

eqtls <- eqtls %>%
  filter(., phenotype_id == "ENSG00000083444.16") %>%
  filter(., variant_id %in% list_all_snps) %>%
  left_join(., rsids, by = "variant_id")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Format tables to use in MRbase analysis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# The effect allele is the second allele in the variant ID
# Need to find a way to systematically include rsid instead of variant_ID

eqtls <- eqtls %>%
  mutate(Phenotype = "PLOD1_eQTL_Whole_Blood",
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

write.table(eqtls, here(project_dir,"MR_analysis","MR_exposures","Whole_Blood_EUR_PLOD1_exposure.tsv"),
            sep = "\t", row.names = F, quote = F)



