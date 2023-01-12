# https://github.com/StoreyLab/qvalue
# script to compute FDR from the multixcan output pvalues
# saves new files in the /home/kidney_genetics_project/multixcan_output/ folder with the prefix "qvalues_"
# can run this script using salloc, as it doesnt require a lot of time/memory
library(qvalue)
library(dplyr)

### BELUGA:

for (QTL in c("eqtl", "sqtl")) {

  home_dir = "/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project/multixcan_output"

  tbl <- read.table(paste0(home_dir,"hematuria_",QTL,"_smultixcan.txt"), header=T, sep= "\t"))

  # compute q-values:
  qvalues <- qvalue(p = tbl$pvalue))
  # remove non-wanted elements:
  qvalues2 <- qvalues[-c(1,2,6,7,8)]

  # merge q-values with original list:
  full_tbl <- cbind(tbl, qvalues2)
  # NOTE: I did not export these results tables with the qvalues, but it is possible to integrate in the for loop below if needed.
  # instead I only output the pvalue threshold that correspond to an FDR < 0.05 in a table, for eqtl and sqtl.

  fdr_threshold <- as.data.frame(matrix(0,1,2))
  colnames(fdr_threshold) <- c("eqtl","sqtl")

  if (QTL=="eqtl") {j = 1}  # condition to save FDR in either eqtl (==1) or sqtl (==2) column 
  else {j = 2}
  full_tbl %>% select("qvalues","pvalues") %>% filter(., qvalues <= 0.01) %>% 
	select("pvalues") %>% max(.) -> fdr_threshold[1,j]  # get FDR < 0.01 threshold and output corresponding pvalue
}

write.table(fdr_threshold, "/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project/multixcan_output/FDR_threshold.tsv", row.names = F, quote = F, sep = "\t")
