# https://github.com/StoreyLab/qvalue
# script to compute FDR from the sprediXcan output pvalues
# saves new files in the /home/renal_genetics_project/spredixcan_output/*qtl/ folder with the prefix "qvalues_"
# can run this script using salloc, as it doesnt require a lot of time/memory
library(qvalue)
library(dplyr)

### BELUGA:

# read list of tissues to output all FDR thresholds:
fdr_threshold <- read.table("/scratch/fridald4/hematuria_scripts/list_of_tissues.txt", header = F)

for (QTL in c("eqtl", "sqtl")) {

  home_dir = paste("/home/fridald4/projects/def-gsarah/fridald4/renal_genetics_project/spredixcan_output/", QTL, sep="")

  # list files:
  filelist_full <- list.files(pattern = "^hematuria_TOPMed_imputed_", path = home_dir, full.names=TRUE)
  # read actual files and save in list
  tbl <- lapply(filelist_full, function(x) read.table(x, header=T, sep= "\t"))
  # remove columns with NA:
  tbl2 <- lapply(tbl, function(x) x[-c(7,8,9)])

  # compute q-values:
  qvalues <- lapply(tbl2, function(x) qvalue(p = x$pvalue))
  # remove non-wanted elements:
  qvalues2 <- lapply(qvalues, function(x) x[-c(1,2,6,7,8)])

  # merge q-values with original list:
  full_tbl <- Map(c, tbl2, qvalues2)
  # NOTE: I did not export these results tables with the qvalues, but it is possible to integrate in the for loop below if needed.
  # instead I only output the pvalue threshold that correspond to an FDR < 0.05 in a table, for eqtl and sqtl.

  for (i in 1:49) {
    if (QTL=="eqtl") {j = 2}  # condition to save FDR in either eqtl (==2) or sqtl (==3) column 
    else {j = 3}
    as.data.frame(full_tbl[[i]]) %>% select("qvalues","pvalues") %>% filter(., qvalues <= 0.01) %>% 
	  select("pvalues") %>% max(.) -> fdr_threshold[i,j]  # get FDR < 0.01 threshold and output corresponding pvalue
  }
}

colnames(fdr_threshold) <- c("tissue", "eqtl_FDR", "sqtl_FDR")
write.table(fdr_threshold, "/home/fridald4/projects/def-gsarah/fridald4/renal_genetics_project/spredixcan_output/FDR_threshold_per_tissue.csv", row.names = F, quote = F, sep = "\t")
