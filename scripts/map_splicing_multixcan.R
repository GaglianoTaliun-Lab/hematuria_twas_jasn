# Map gene names from multixcan output to chr and position.

library(dplyr)
library(tidyr)

HOMEPATH="/home/fridald4/projects/def-gsarah/fridald4/hematuria_project/multixcan_output/"

#### for sqtl, there are no gene names, therefore need to separate the chr and start-end position from gene_name

sqtl_results = read.table(paste0(HOMEPATH,"hematuria_sqtl_smultixcan.txt"), sep = "\t", header = TRUE)
sqtl_results$gene_name <- gsub("intron_", "", sqtl_results$gene_name)
sqtl_results %>%
  separate(gene_name, c("chromosome", "position_start", "position_end")) -> sqtl_results2

### write output into new file:
out_dir = "/home/fridald4/projects/def-gsarah/fridald4/hematuria_project/multixcan_output/"
write.table(sqtl_results2, paste0(out_dir, "hematuria_sqtl_smultixcan_mapped.txt"), sep = "\t", row.names = F, quote = F)

