# Description: run colocalization analysis across two traits

# Packages -------------------------------------------------------

library(coloc)
library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(colochelpR)

# Arguments ----------------------------------------------------------------------

project_dir <- "/home/fridald4/projects/def-gsarah/fridald4/renal_genetics_project/MR_analysis"

p1 = 1e-04
p2 = 1e-04
p12 = 1e-05

# Read files -------------------------------------------------------

input_coloc_gwas <- read.table(here(project_dir, "colocalization", "data_for_coloc", "hematuria_ukbb_P3H2.tsv"), sep = "\t", header  = T)
input_coloc_qtls <- read.table(here(project_dir, "colocalization", "data_for_coloc", "P3H2_qtls.tsv"), sep = "\t", header = T)

# Run coloc ----------------------------------------------------------------------------------

df2 <- input_coloc_qtls %>%
  filter(., maf > 0) %>%
  rename(MAF = maf, pvalues = p.value) 

df1 <- input_coloc_gwas %>%
  filter(., maf > 0) %>%
  filter(!duplicated(SNP)) %>%
  rename(MAF = maf, pvalues = p.value)
    
coloc_results <- coloc.abf(dataset1 = list(type = "quant",
                                           snp = df1$SNP,
                                           beta = df1$beta,
                                           varbeta = df1$varbeta,
					   pvalues = df1$pvalues,
                                           MAF = df1$MAF,
                                           N = df1$N),
                           dataset2 = list(type = "quant",
                                           snp = df2$SNP,
                                           beta = df2$beta,
                                           varbeta = df2$varbeta,
					   pvalues = df2$pvalues,
                                           MAF = df2$MAF,
                                           N = df2$N),
                           p1 = p1, p2 = p2, p12 = p12)
      
coloc_results_summ <- coloc_results$summary
coloc_results_res <- coloc_results$results
      
write.table(coloc_results_summ, here(project_dir, "colocalization", str_c("coloc_summary_", df2$eQTL_dataset[1],"_",df1$GWAS[1],".tsv")), sep = "\t", row.names = T, quote = F, col.names = F)
write.table(coloc_results_res, here(project_dir, "colocalization", str_c("coloc_results_", df2$eQTL_dataset[1],"_",df1$GWAS[1],".tsv")), sep = "\t", row.names = F, quote = F) 
