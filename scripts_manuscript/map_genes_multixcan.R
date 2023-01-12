# Map gene names from multixcan output to chr and position.

library(dplyr)

#### for sqtl, there are no gene names, therefore need to separate the chr and start-end position from gene_name

home_dir = "/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project/multixcan_output/"

### read reference from Homo sapiens obtained from: https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
### this file includes gene symbols and IDs and chromosome, but not start/end position of the gene
read.table("/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project/gene_info_chr_maplocation.txt", sep = "\t", header = T) %>%
  filter(., chromosome %in% (1:22)) -> gene_symbol_reference

### read reference from Homo sapiens obtained from: https://ftp.ncbi.nih.gov/gene/DATA/gene_neighbors.gz
### this file includes gene IDs, chromosome, start/end positions, but not gene symbols.
read.table("/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project/geneID_pos_maplocation.txt", sep = "\t", header = T) %>%
  filter(., chromosome %in% (1:22)) -> gene_ID_reference

### read output from multixcan:
multixcan <- read.table(paste0(home_dir,"hematuria_eqtl_smultixcan.txt"), sep = "\t", header = T)

### join gene_symbol table with multixcan output to include chr and geneID:
genes_mapped <- left_join (multixcan, gene_symbol_reference, by = c("gene_name" = "Symbol"))

### join gene_ID table with genes_mapped list to include start/end positions of gene:
genes_mapped2 <- left_join (genes_mapped, gene_ID_reference, by = "GeneID")

genes_mapped2 %>%
  dplyr::mutate(chromosome = chromosome.x) %>%
  dplyr::select(
    gene,
    gene_name,
    pvalue,
    n,
    n_indep,
    p_i_best,
    t_i_best,
    p_i_worst,
    t_i_worst,
    eigen_max,
    eigen_min,
    eigen_min_kept,
    z_min,
    z_max,
    z_mean,
    z_sd,
    tmi,
    status,
    GeneID,
    chromosome,
    map_location,
    start_position,
    end_position) -> genes_mapped3

### write output into new file:
out_dir = "/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project/multixcan_output/"
write.table(genes_mapped3, paste0(out_dir, "hematuria_eqtl_smultixcan_mapped.txt"), sep = "\t", row.names = F, quote = F)

