# Map gene names from spredixcan output to chr and position.

library(dplyr)

#### for sqtl, there are no gene names, therefore need to separate the chr and start-end position from gene_name

home_dir = paste("/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project/spredixcan_output/eqtl/",sep="")

### read reference from Homo sapiens obtained from: https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
### this file includes gene symbols and IDs and chromosome, but not start/end position of the gene
read.table("/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project/gene_info_chr_maplocation.txt", sep = "\t", header = T) %>%
  filter(., chromosome %in% (1:22)) -> gene_symbol_reference

### read reference from Homo sapiends obtained from: https://ftp.ncbi.nih.gov/gene/DATA/gene_neighbors.gz
### this file includes gene IDs, chromosome, start/end positions, but not gene symbols.
read.table("/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project/geneID_pos_maplocation.txt", sep = "\t", header = T) %>%
  filter(., chromosome %in% (1:22)) -> gene_ID_reference

### read outputs from predixcan and save in list:
filelist_full <- list.files(pattern = "^hematuria_TOPMed_imputed_eqtl_", path = home_dir, full.names=TRUE)
genes_per_tissue <- lapply(filelist_full, function(x) read.table(x, header= TRUE))

### join gene_symbol table with spredixcan output to include chr and geneID:
genes_mapped <- lapply(genes_per_tissue, function(x)left_join (x, gene_symbol_reference, by = c("gene_name" = "Symbol")))

### join gene_ID table with genes_mapped list to include start/end positions of gene:
genes_mapped2 <- lapply(genes_mapped, function(x)left_join (x, gene_ID_reference, by = "GeneID"))

### write output into new files - one dataframe per tissue:
out_dir = "/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project/spredixcan_output/gene_mapped_results/"
filelist_out <- list.files(pattern = "^hematuria_TOPMed_imputed_eqtl_", path = home_dir, full.names=FALSE)

for (i in 1:49) {
  as.data.frame(genes_mapped2[[i]]) %>%
    .[,-c(7, 8, 9, 15, 18)] %>%  # remove NA columns from spredixcan output
    write.table(., paste(out_dir,"mapped_", filelist_out[i],sep=""), sep = "\t", row.names = F, quote = F)
}
