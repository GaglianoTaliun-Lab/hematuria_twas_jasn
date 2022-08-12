### script to count total number of genes in MASHR models (eqtl and sqtl) in order to obtain exact proportion
### of genes/splicing sites that are overlapped and were used in the TWAS per tissue.

library(RSQLite)
library(dplyr)
library(ggplot2)

# set home directory
home_dir = "/home/fridald4/projects/def-gsarah/fridald4/hematuria_project/"

# read tissues list:
tissues_list <- read.table("/scratch/fridald4/hematuria_scripts/list_of_tissues.txt")
colnames(tissues_list) <- "tissue"

# create empty output data frames:
number_genes <- as.data.frame(matrix(data=0, nrow=49, ncol=1))
number_sites <- as.data.frame(matrix(data=0, nrow=49, ncol=1))

# loop start:
for (i in 1:49) {
  filename <- tissues_list[i,1]
  sqlite.driver <- dbDriver("SQLite")
  
  # for eqtl:
  db <- dbConnect(sqlite.driver,
                  dbname = paste(home_dir,"MASHR/eqtl/mashr/mashr_", filename, ".db", sep=""))
                  dbListTables(db)
                  extra_info_genes <- dbReadTable(db,"extra")
                  number_genes[i,1] <- nrow(extra_info_genes)
  # for sqtl:
  db <- dbConnect(sqlite.driver,
                  dbname = paste(home_dir,"MASHR/sqtl/mashr/mashr_", filename,".db", sep=""))
                  extra_info_sites <- dbReadTable(db,"extra")
                  number_sites[i,1] <- nrow(extra_info_sites)
}

# dbDisconnect()

# bind tissues names and number of genes/sites:
colnames(number_genes) <- "total_genes"
colnames(number_sites) <- "total_sites"

number_genes_per_tissue <- cbind(tissues_list, number_genes)
number_sites_per_tissue <- cbind(tissues_list, number_sites)

write.table(number_genes_per_tissue, paste(home_dir,"MASHR/number_total_genes_per_tissue.txt",sep=""), quote = F, sep = " ",row.names = F)
write.table(number_sites_per_tissue, paste(home_dir,"MASHR/number_total_sites_per_tissue.txt",sep=""), quote = F, sep = " ",row.names = F)

###### I have obtained the total number of genes/splicing sites that were output in the TWAS results by counting the number of rows.
###### To obtain the proportion of those, compared to the total number of genes/sites available in the models:

# read csv file:
twas_genes_sites <- read.csv(paste(home_dir, "spredixcan_output/summary_results_spredixcan_hematuria.csv", sep=""), header=T)[,1:5]

# join tables:
full_join(number_genes_per_tissue, number_sites_per_tissue, by = "tissue") %>%
  full_join(., twas_genes_sites) -> proportion_genes

# get %:
proportion_genes["proportion_genes_in_twas"] = (proportion_genes$eQTL_number_of_genes / proportion_genes$total_genes) * 100
proportion_genes["proportion_sites_in_twas"] = (proportion_genes$sQTL_number_of_splicing_sites / proportion_genes$total_sites) * 100

# save table:
write.table(proportion_genes, paste(home_dir,"spredixcan_output/figures/proportion_of_genes_and_sites_per_tissue_49tissues.txt",sep=""), quote = F, sep = "\t",row.names = F)

# save table after removing sex-specific tissues:
filter(proportion_genes, tissue != "Ovary") %>%
  filter(., tissue != "Testis") %>%
  filter(., tissue != "Uterus") %>%
  filter(., tissue != "Prostate") %>%
  filter(., tissue != "Vagina") -> proportion_genes_44tissues

# save table:
write.table(proportion_genes_44tissues, paste(home_dir,"spredixcan_output/figures/proportion_of_genes_and_sites_per_tissue_44tissues.txt",sep=""), quote = F, sep = "\t",row.names = F)

#### plot proportion of genes:
ggplot(data = proportion_genes, mapping = aes(x=tissue, y=proportion_genes_in_twas)) +
  geom_col() + 
  theme_bw() +
  # guides(x = guide_axis(n.dodge = 2), y.sec = guide_axis())
  guides(x = guide_axis(angle = 90))
ggsave(paste(home_dir, "spredixcan_output/figures/proportion_of_genes.pdf", sep=""))

#### plot proportion of splicing sites:
ggplot(data = proportion_genes, mapping = aes(x=tissue, y=proportion_sites_in_twas)) +
  geom_col() + 
  theme_bw() +
  # guides(x = guide_axis(n.dodge = 2), y.sec = guide_axis())
  guides(x = guide_axis(angle = 90))
ggsave(paste(home_dir, "spredixcan_output/figures/proportion_of_sites.pdf", sep=""))

###### Finally, create the same plots of proportions above but first remove the sex-specific tissues:
###### ovary, testis, uterus, prostate, vagina

filter(proportion_genes, tissue != "Ovary") %>%
  filter(., tissue != "Testis") %>%
  filter(., tissue != "Uterus") %>%
  filter(., tissue != "Prostate") %>%
  filter(., tissue != "Vagina") -> proportion_genes_44

# write table:
write.table(proportion_genes_44, paste(home_dir,"spredixcan_output/figures/proportion_of_genes_and_sites_per_tissue_44tissues.txt",sep=""), quote = F, sep = "\t",row.names = F)

#### plot proportion of genes:
ggplot(data = proportion_genes_44, mapping = aes(x=tissue, y=proportion_genes_in_twas)) +
  geom_col() + 
  theme_bw() +
  # guides(x = guide_axis(n.dodge = 2), y.sec = guide_axis())
  guides(x = guide_axis(angle = 90))
ggsave(paste(home_dir, "spredixcan_output/figures/proportion_of_genes.pdf", sep=""))

#### plot proportion of splicing sites:
ggplot(data = proportion_genes_44, mapping = aes(x=tissue, y=proportion_sites_in_twas)) +
  geom_col() + 
  theme_bw() +
  # guides(x = guide_axis(n.dodge = 2), y.sec = guide_axis())
  guides(x = guide_axis(angle = 90))
ggsave(paste(home_dir, "spredixcan_output/figures/proportion_of_sites.pdf", sep=""))

### import file with sample size per tissue:
tissues_sample_size <- read.csv(paste(home_dir,"GTEx_number_samples.csv",sep=""), header = T)

# merge with proportion_genes_44 and arrange by sample size:
inner_join(proportion_genes_44, tissues_sample_size, by = "tissue") %>%
  arrange(., samples) -> tissues_sample_size_44

### number of significant genes per tissue
ggplot(data = proportion_genes_44, mapping = aes(x=tissue, y=significant_genes)) +
  geom_col() + 
  theme_bw() +
  # guides(x = guide_axis(n.dodge = 2), y.sec = guide_axis())
  guides(x = guide_axis(angle = 90))
ggsave(paste(home_dir, "spredixcan_output/figures/eQTL_significant_genes.pdf", sep=""))

### number of significant splicing sites per tissue
ggplot(data = proportion_genes_44, mapping = aes(x=tissue, y=significant_splicing_sites)) +
  geom_col() + 
  theme_bw() +
  # guides(x = guide_axis(n.dodge = 2), y.sec = guide_axis())
  guides(x = guide_axis(angle = 90))
ggsave(paste(home_dir, "spredixcan_output/figures/sQTL_significant_splicing_sites.pdf", sep=""))
