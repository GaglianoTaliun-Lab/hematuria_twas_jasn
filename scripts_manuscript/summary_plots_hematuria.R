## Script to create summary plots of general TWAS results.

## The first section plots all 49 tissues from the GTEx, ordered by increasing number of sample sizes.
## The second section plots 44 tissues (removing sex-specific tissues: ovary, testis, prostate, vagina and uterus),
## also in increasing number of sample sizes.
## On the third section I calculate the correlations between genes/sites and the sample size and output a table with results.

## script ran on local computer, but it only uses one file as input (see below), therefore it is easily transferrable. 

# Import libraries and required files:

library(dplyr)
library(ggplot2)

setwd("~/Documents/research-projects/hematuria_project")

summary <- read.delim("proportion_of_genes_and_sites_per_tissue.txt", header=T)
## note: this file was created on Béluga cluster by merging the csv file in local computer (summary_results_spredixcan_hematuria.csv)
## and the stats obtained through the "count_number_of_genes_mash.R" script.

### import file with sample size per tissue:
tissues_sample_size <- read.csv("GTEx_number_samples.csv", header = T)

####### merge with summary data frame and arrange by sample size:
inner_join(summary, tissues_sample_size, by = "tissue") %>%
  arrange(., samples) -> summary_arranged


######################################## SECTION 1

# number of total genes per tissue
ggplot(data = summary_arranged, mapping = aes(x= reorder(tissue, -samples), y=eQTL_number_of_genes)) +
  geom_col() + 
  theme_bw() +
  guides(x = guide_axis(angle = 90)) +
  labs(x = "GTEx Tissue", y = "Number of total genes in TWAS")
ggsave("eQTL_number_of_genes.pdf")

# number of total splicing sites per tissue
ggplot(data = summary_arranged, mapping = aes(x=reorder(tissue, -samples), y=sQTL_number_of_splicing_sites)) +
  geom_col() + 
  theme_bw() +
  guides(x = guide_axis(angle = 90)) +
  labs(x = "GTEx Tissue", y = "Number of total splicing sites in TWAS")
ggsave("sQTL_number_of_splicing_sites.pdf")

# number of significant genes per tissue
ggplot(data = summary_arranged, mapping = aes(x=reorder(tissue, -samples), y=significant_genes)) +
  geom_col() + 
  theme_bw() +
  guides(x = guide_axis(angle = 90)) +
  labs(x = "GTEx Tissue", y = "Number of significant genes in TWAS")
ggsave("eQTL_significant_genes_49tissues.pdf")

# number of significant splicing sites per tissue
ggplot(data = summary_arranged, mapping = aes(x=reorder(tissue, -samples), y=significant_splicing_sites)) +
  geom_col() + 
  theme_bw() +
  guides(x = guide_axis(angle = 90)) +
  labs(x = "GTEx Tissue", y = "Number of significant splicing sites in TWAS")
ggsave("sQTL_significant_splicing_sites_49tissues.pdf")

## proportion of genes compared to total genes in model:
ggplot(data = summary_arranged, mapping = aes(x=reorder(tissue, -samples), y=proportion_genes_in_twas)) +
  geom_col() + 
  theme_bw() +
  guides(x = guide_axis(angle = 90)) +
  labs(x = "GTEx Tissue", y = "Percentage of TWAS genes")
ggsave("proportion_of_genes_49tissues.pdf")

## proportion of splicing sites compared to total sites in model:
ggplot(data = summary_arranged, mapping = aes(x=reorder(tissue, -samples), y=proportion_sites_in_twas)) +
  geom_col() + 
  theme_bw() +
  guides(x = guide_axis(angle = 90)) +
  labs(x = "GTEx Tissue", y = "Percentage of TWAS splicing sites")
ggsave("proportion_of_sites_49tissues.pdf")


########################################## SECTION 2

###### Same plots as above but first remove the sex-specific tissues:
###### ovary, testis, uterus, prostate, vagina

filter(summary_arranged, tissue != "Ovary") %>%
  filter(., tissue != "Testis") %>%
  filter(., tissue != "Uterus") %>%
  filter(., tissue != "Prostate") %>%
  filter(., tissue != "Vagina") -> summary_44_arranged

# number of total genes per tissue
ggplot(data = summary_44_arranged, mapping = aes(x=reorder(tissue, -samples), y=eQTL_number_of_genes)) +
  geom_col() + 
  theme_bw() +
  guides(x = guide_axis(angle = 90)) +
  labs(x = "GTEx Tissue", y = "Number of total genes in TWAS")
ggsave("eQTL_number_of_genes_44tissues.pdf")

# number of total splicing sites per tissue
ggplot(data = summary_44_arranged, mapping = aes(x=reorder(tissue, -samples), y=sQTL_number_of_splicing_sites)) +
  geom_col() + 
  theme_bw() +
  guides(x = guide_axis(angle = 90)) +
  labs(x = "GTEx Tissue", y = "Number of total splicing sites in TWAS")
ggsave("sQTL_number_of_splicing_sites_44tissues.pdf")

# number of significant genes per tissue
ggplot(data = summary_44_arranged, mapping = aes(x=reorder(tissue, -samples), y=significant_genes)) +
  geom_col() + 
  theme_bw() +
  guides(x = guide_axis(angle = 90)) +
  labs(x = "GTEx Tissue", y = "Number of significant genes in TWAS")
ggsave("eQTL_significant_genes_44tissues.pdf")

# number of significant splicing sites per tissue
ggplot(data = summary_44_arranged, mapping = aes(x=reorder(tissue, -samples), y=significant_splicing_sites)) +
  geom_col() + 
  theme_bw() +
  guides(x = guide_axis(angle = 90)) +
  labs(x = "GTEx Tissue", y = "Number of significant splicing sites in TWAS")
ggsave("sQTL_significant_splicing_sites_44tissues.pdf")

#### plot proportion of genes compared to all genes in model:
ggplot(data = summary_44_arranged, mapping = aes(x=reorder(tissue, -samples), y=proportion_genes_in_twas)) +
  geom_col() + 
  theme_bw() +
  guides(x = guide_axis(angle = 90)) +
  labs(x = "GTEx Tissue", y = "Percentage of TWAS genes")
ggsave("proportion_of_genes_44tissues.pdf")

#### plot proportion of splicing sites compared to all sites in model:
ggplot(data = summary_44_arranged, mapping = aes(x=reorder(tissue, -samples), y=proportion_sites_in_twas)) +
  geom_col() + 
  theme_bw() +
  guides(x = guide_axis(angle = 90)) +
  labs(x = "GTEx Tissue", y = "Percentage of TWAS splicing sites")
ggsave("proportion_of_sites_44tissues.pdf")

############################################ SECTION 3

### correlations using 49 tissues:

correlations_49tissues <- data.frame(matrix(0,8,2), 
                                     row.names = c("total_genes", "total_sites", "TWAS_number_of_genes",
                                                   "TWAS_number_of_sites", "significant_genes", "significant_sites",
                                                   "proportion_genes_TWAS", "proportion_sites_TWAS"))
colnames(correlations_49tissues) <- c("cor_value", "pvalue")

for (i in 2:9) {
  k = i-1
  correlations_49tissues[k,1] <- cor.test(summary_arranged[,i], summary_arranged$samples) %>% .$estimate %>% round(., 3)
  correlations_49tissues[k,2] <- cor.test(summary_arranged[,i], summary_arranged$samples) %>% .$p.value %>% format(., scientific = T)
}

write.table(correlations_49tissues, "correlations_GTEx_sample_size_and_TWAS_genes_49tissues.csv", quote = F, sep = ",", row.names = T, col.names = NA)

## correlations using 44 tissues:

correlations_44tissues <- data.frame(matrix(0,8,2), 
                                     row.names = c("total_genes", "total_sites", "TWAS_number_of_genes",
                                                   "TWAS_number_of_sites", "significant_genes", "significant_sites",
                                                   "proportion_genes_TWAS", "proportion_sites_TWAS"))
colnames(correlations_44tissues) <- c("cor_value", "pvalue")

for (i in 2:9) {
  k = i-1
  correlations_44tissues[k,1] <- cor.test(summary_44_arranged[,i], summary_44_arranged$samples) %>% .$estimate %>% round(., 3)
  correlations_44tissues[k,2] <- cor.test(summary_44_arranged[,i], summary_44_arranged$samples) %>% .$p.value %>% format(., scientific = T)
}

write.table(correlations_44tissues, "correlations_GTEx_sample_size_and_TWAS_genes_44tissues.csv", quote = F, sep = ",", row.names = T, col.names = NA)