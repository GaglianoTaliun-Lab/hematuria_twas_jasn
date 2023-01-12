# Description - obtain outcomes for the PLOD2 gene across different GWAS summary statistics
                                               
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load libraries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(R.utils)
library(tidyverse)
library(data.table)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir <- here("/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project")

snps_main <- c("rs7641322","rs57655223","rs73148997")
snps_proxy_palindromic <- c("rs1170389")
snps_proxy_rs7641322 <- c("rs1449444")
snps_proxy_rs73148997 <- c("rs4681295")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read outcome files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Hematuria UKBB (TOPMed imputed) GWAS summary statistics:
hematuria_UKBB <- fread(here(project_dir,"summary_statistics","hematuria_sumstats_with_chrpos_ID.txt.gz"))

# ACR UKBB GWAS summary statistics:
ACR_UKBB <- fread(here(project_dir,"MR_analysis","ACR_sumstats","ratio.ukb_chr3_v3.SAIGE.txt"))

# Composite phenotype (hematuria + ACR form UKBB):
composite_UKBB <- fread(here(project_dir,"MR_analysis","composite","ratio.casecontrol.ukb_chr3_v3.SAIGE.txt"))

# Stanzick et al 2021 eGFR cystatin meta-analysis between UKBB and CKDGEN:
meta_all_eGFRcrea_stanzick2021 <- fread(here(project_dir,"MR_analysis","stanzick2021_meta_analysis_eGFRcrea_CKDGEN_UKBB.gc"))

# Stanzick et al 2021 EUR eGFR creatinine meta-analysis between UKBB and CKDGEN:
meta_eGFRcrea_stanzick2021 <- fread(here(project_dir,"MR_analysis","stanzick2021_meta_analysis_EUR_eGFRcrea_CKDGEN_UKBB.gc"))

# Stanzick et al 2021 EUR eGFR cystatin meta-analysis between UKBB and CKDGEN:
meta_eGFRcys_stanzick2021 <- fread(here(project_dir,"MR_analysis","stanzick2021_meta_analysis_EUR_eGFRcys_CKDGEN_UKBB.gc"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract SNPs of interest
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

list_all_snps <- c(snps_main, snps_proxy_palindromic, snps_proxy_rs7641322, snps_proxy_rs73148997)

hematuria_UKBB <- hematuria_UKBB %>%
  filter(., rsids %in% list_all_snps)

ACR_UKBB <- ACR_UKBB %>%
  filter(., rsid %in% list_all_snps)

composite_UKBB <- composite_UKBB %>%
  filter(., rsid %in% list_all_snps)

meta_all_eGFRcrea_stanzick2021 <- meta_all_eGFRcrea_stanzick2021 %>%
  filter(., RSID %in% list_all_snps)

meta_eGFRcrea_stanzick2021 <- meta_eGFRcrea_stanzick2021 %>%
  filter(., RSID %in% list_all_snps)

meta_eGFRcys_stanzick2021 <- meta_eGFRcys_stanzick2021 %>%
  filter(., RSID %in% list_all_snps)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Format tables to use in MRbase analysis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

hematuria_UKBB <- hematuria_UKBB %>%
  mutate(N = 394591,
	 phenotype = "hematuria_UKBB_phecode593") %>%
  select(
    chr = chrom,
    pos,
    other_allele = ref,
    effect_allele = alt,
    SNP = rsids,
    phenotype,
    pval = pval,
    beta = beta,
    se = sebeta,
    EAF = af,
    N
  )

ACR_UKBB <- ACR_UKBB %>%
  mutate(phenotype = "ACR_UKBB") %>%
  select(
    chr = CHR,
    pos = POS,
    other_allele = Allele1,
    effect_allele = Allele2,
    SNP = rsid,
    phenotype,
    pval = p.value,
    beta = BETA,
    se = SE,
    EAF = AF_Allele2,
    N
  )

composite_UKBB <- composite_UKBB %>%
  mutate(N = 345938,
	 phenotype = "composite_UKBB") %>%
  select(
    chr = CHR,
    pos = POS,
    other_allele = Allele1,
    effect_allele = Allele2,
    SNP = rsid,
    phenotype,
    pval = p.value,
    beta = BETA,
    se = SE,
    EAF = AF_Allele2,
    N
  )

meta_all_eGFRcrea_stanzick2021 <- meta_all_eGFRcrea_stanzick2021 %>%
  mutate(N = 1201929,
         phenotype = "meta_all_eGFRcrea_stanzick2021") %>%
  select(
    chr,
    pos,
    other_allele = Allele2,
    effect_allele = Allele1,
    SNP = RSID,
    phenotype,
    pval = P.value.GC,
    beta = Effect,
    se = StdErr.GC,
    EAF = Freq1,
    N
  )


meta_eGFRcrea_stanzick2021 <- meta_eGFRcrea_stanzick2021 %>%
  mutate(N = 1004040,
	 phenotype = "meta_EUR_eGFRcrea_stanzick2021") %>%
  select(
    chr,
    pos,
    other_allele = Allele2,
    effect_allele = Allele1,
    SNP = RSID,
    phenotype,
    pval = P.value.GC,
    beta = Effect,
    se = StdErr.GC,
    EAF = Freq1,
    N
  )

meta_eGFRcys_stanzick2021 <- meta_eGFRcys_stanzick2021 %>%
  mutate(N = 460826,
         phenotype = "meta_EUR_eGFRcys_stanzick2021") %>%
  select(
    chr,
    pos,
    other_allele = Allele2,
    effect_allele = Allele1,
    SNP = RSID,
    phenotype,
    pval = P.value.GC,
    beta = Effect,
    se = StdErr.GC,
    EAF = Freq1,
    N
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write outcome files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.table(hematuria_UKBB, here(project_dir,"MR_analysis","MR_outcomes","hematuria_TOPMed_imputed_PLOD2_outcome.tsv"),
            sep = "\t", row.names = F, quote = F)

write.table(ACR_UKBB, here(project_dir,"MR_analysis","MR_outcomes","ACR_UKBB_PLOD2_outcome.tsv"),
            sep = "\t", row.names = F, quote = F)

write.table(composite_UKBB, here(project_dir,"MR_analysis","MR_outcomes","composite_hematuria_ACR_UKBB_PLOD2_outcome.tsv"),
            sep = "\t", row.names = F, quote = F)

write.table(meta_all_eGFRcrea_stanzick2021, here(project_dir,"MR_analysis","MR_outcomes","metaanalysis_stanzick2021_eGFRcreatinine_PLOD2_outcome.tsv"),
            sep = "\t", row.names = F, quote = F)

write.table(meta_eGFRcrea_stanzick2021, here(project_dir,"MR_analysis","MR_outcomes","metaanalysis_stanzick2021_EUR_eGFRcreatinine_PLOD2_outcome.tsv"),
            sep = "\t", row.names = F, quote = F)

write.table(meta_eGFRcys_stanzick2021, here(project_dir,"MR_analysis","MR_outcomes","metaanalysis_stanzick2021_EUR_eGFRcystatin_PLOD2_outcome.tsv"),
            sep = "\t", row.names = F, quote = F)

