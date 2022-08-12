# Description: run colocalization analysis across kidney disease traits on genes PLOD1 and PLOD2

# 1. PLOD1 expression and eGFRcrea (meta-analysis)
# 2. PLOD1 expression and eGFRcys (meta-analysis)
# 3. PLOD2 expression and eGFRcys (CKDGen)

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

GWAS_for_coloc <- list.files(here(project_dir, "colocalization", "data_for_coloc"), pattern = "^eGFR", all.files = T, full.names = T) %>% sort()
GWASnames_for_coloc <- list.files(here(project_dir, "colocalization", "data_for_coloc"), pattern = "^eGFR", all.files = T, full.names = F) %>% sort() %>%
  str_remove(., ".tsv")

QTLs_for_coloc <- list.files(here(project_dir, "colocalization", "data_for_coloc"), pattern = "*qtls.tsv", all.files = T, full.names = T) %>% sort()
QTLnames_for_coloc <- list.files(here(project_dir, "colocalization", "data_for_coloc"), pattern = "*qtls.tsv", all.files = T, full.names = F) %>% sort() %>%
  str_remove(., ".tsv")

input_coloc_gwas <- list()
input_coloc_qtls <- list()

input_coloc_gwas <- setNames(GWAS_for_coloc %>%
                           lapply(., function(x) read.table(x, sep = "\t", header  = T)),
                         nm = GWASnames_for_coloc
)

input_coloc_qtls <- setNames(QTLs_for_coloc %>%
                               lapply(., function(x) read.table(x, sep = "\t", header  = T)),
                             nm = QTLnames_for_coloc
)

# Run coloc ----------------------------------------------------------------------------------

coloc_results_summ <- list()
coloc_results_res <- list()
results_names <- array()
i=1

for (qtl in 1:length(QTLs_for_coloc)) {
  for (gwas in 1:length(GWAS_for_coloc)) {
    
    df2 <- input_coloc_qtls[[qtl]] %>%
      filter(., maf > 0) %>%
      rename(MAF = maf, pvalues = p.value)
    df1 <- input_coloc_gwas[[gwas]] %>%
      filter(., maf > 0) %>%
      rename(MAF = maf, pvalues = p.value)
    
    if (str_detect(GWAS_for_coloc[gwas], "PLOD1") == TRUE && str_detect(QTLs_for_coloc[qtl], "PLOD1") == TRUE) {
      
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
                                    p1 = p1, p2 = p2, p12 = p12
      )
      
      coloc_results_summ[[i]] <- coloc_results$summary
      coloc_results_res[[i]] <- coloc_results$results
      results_names[i] <- str_c(df1$GWAS[1],"_",df2$eQTL_dataset[1])
      i=i+1
      
    }
    
    else if (str_detect(GWAS_for_coloc[gwas], "PLOD2") == TRUE && str_detect(QTLs_for_coloc[qtl], "PLOD2") == TRUE) {
      
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
                                    p1 = p1, p2 = p2, p12 = p12
      )
      
      coloc_results_summ[[i]] <- coloc_results$summary
      coloc_results_res[[i]] <- coloc_results$results
      results_names[i] <- str_c(df1$GWAS[1],"_",df2$eQTL_dataset[1])
      i=i+1
      
    } else next
  
  }
  
}

for (i in 1:length(results_names)) {

  summ_out <- as.data.frame(coloc_results_summ[[i]])
  res_out <- as.data.frame(coloc_results_res[[i]])
  write.table(summ_out, here(project_dir, "colocalization", str_c("coloc_summary_", results_names[i],".tsv")), sep = "\t", row.names = T, quote = F, col.names = F)
  write.table(res_out, here(project_dir, "colocalization", str_c("coloc_results_", results_names[i],".tsv")), sep = "\t", row.names = F, quote = F) 

}
