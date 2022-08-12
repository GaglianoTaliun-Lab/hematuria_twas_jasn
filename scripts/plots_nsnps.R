library(dplyr)
library(ggplot2)

setwd("~/Documents/research-projects/hematuria_project")

# choose to plot either eQTL or sQTL by defining QTL (the rest is automated):
QTL="eQTL"
QTL="sQTL"

# choose source of input data:
source="predixcan"
#source="dbtables"

tissue_number <- read.table("GTEx_number_samples.csv", sep = ",", header = T)

if (QTL == "eQTL" && source == "dbtables") {
  tissue_df <- read.table("snps_in_models/snps_in_genes_per_tissue.tsv", sep = "\t", header = T)
  type = "gene"
} else if (QTL == "sQTL" && source == "dbtables") {
  tissue_df <- read.table("snps_in_models/snps_in_sites_per_tissue.tsv", sep = "\t", header = T)
  type="site" 
} else if (QTL == "eQTL" && source == "predixcan") {
  read.table("snps_in_models/snps_in_genes_per_tissue_from_spredixcan_results.tsv", sep = "\t", header = T) %>%
    left_join(., tissue_number, by = "tissue") -> tissue_df
  colnames(tissue_df) <- c("tissue","gene","gene_name","n_snps","samples_total","samples_eur")
  type = "gene"
} else if (QTL == "sQTL" && source == "predixcan") {
  tissue_df <- read.table("snps_in_models/snps_in_sites_per_tissue_from_spredixcan_results.tsv", sep = "\t", header = T) %>%
    left_join(., tissue_number, by = "tissue") -> tissue_df
  colnames(tissue_df) <- c("tissue","gene","n_snps","samples_total","samples_eur")
  type="site" 
}

# remove sex-specific tissues:
# filter(tissue_df, tissue != "Ovary") %>%
#   filter(., tissue != "Testis") %>%
#   filter(., tissue != "Uterus") %>%
#   filter(., tissue != "Prostate") %>%
#   filter(., tissue != "Vagina") %>%
#   filter(., n_snps > 0) -> tissue_df

####------------ plot histogram (proportions stratified by nsnps)
ggplot(data = tissue_df, mapping = aes(x= reorder(tissue, -samples_eur), 
                                            fill = as.factor(n_snps))) +
  geom_bar(position = "fill") + 
  theme_bw() +
  guides(x = guide_axis(angle = 90), fill=guide_legend(nrow=1, byrow=TRUE)) +
  scale_fill_manual(values = c("#1B9E77","#D95F02","#7570B3","#E7298A",
                    "#66A61E","#E6AB02","#A6761D","#666666",
                    "red","blue","black","purple"), name = "Number of SNPs") +
  labs(x = "GTEx Tissue", y = paste("Proportion of ",QTL,"s",sep="")) +
  theme(legend.position="bottom", legend.text = element_text(size=8), legend.title = element_text(size=10),
        axis.text.x = element_text(size = 14, angle = 90, hjust = .5, vjust = .5),
        axis.text.y = element_text(size = 14, hjust = .5, vjust = .5),
        axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18))
if (source == "predixcan") {
  ggsave(paste("./figures/summary_results/proportion_",QTL,"s_per_gene_across_tissues_stratified_by_nsnps_predixcan_results.png",sep=""), dpi = 300)
} else if (source == "dbtables") {
  ggsave(paste("./figures/summary_results/proportion_",QTL,"s_per_gene_across_tissues_stratified_by_nsnps_dbtables.png",sep=""), dpi = 300)
}
  
#####------------ plot histogram (total stratified by snps)
ggplot(data = tissue_df, mapping = aes(x= reorder(tissue, -samples_eur), 
                                             fill = as.factor(n_snps))) +
  geom_bar() + 
  theme_bw() +
  guides(x = guide_axis(angle = 90)) +
  scale_fill_manual(values = c("#1B9E77","#D95F02","#7570B3","#E7298A",
                               "#66A61E","#E6AB02","#A6761D","#666666",
                               "red","blue","black","purple"), name = "Number of SNPs") +
  labs(x = "GTEx Tissue", y = paste("Number of ",QTL,"s \n across all genes",sep=""))
if (source == "predixcan") {
ggsave(paste("./figures/summary_results/total_",QTL,"s_per_gene_across_tissues_stratified_by_nsnps_predixcan_results.png",sep=""), dpi = 300)
} else if (source == "dbtables") {
  ggsave(paste("./figures/summary_results/total_",QTL,"s_per_gene_across_tissues_stratified_by_nsnps_dbtables.png",sep=""), dpi = 300)
}
  
#####------------ plot histogram (total without stratification)
ggplot(data = tissue_df, mapping = aes(x= reorder(tissue, -samples_eur))) +
  geom_bar() + 
  theme_bw() +
  guides(x = guide_axis(angle = 90)) +
  geom_text(aes(label = ..count..), stat = "count", vjust = -1, colour = "black", size = 2.5) +
  labs(x = "GTEx Tissue", y = paste("Number of ",QTL,"s \n across all genes",sep=""))
if (source == "predixcan") {
ggsave(paste("./figures/summary_results/total_",QTL,"s_across_tissues_predixcan_results.png",sep=""), dpi = 300)
} else if (source == "dbtables") {
  ggsave(paste("./figures/summary_results/total_",QTL,"s_across_tissues_dbtables.png",sep=""), dpi = 300)
}
  
#####------------ plot histogram with n_snps > 1
tissue_df_nsnp_2 <- filter(tissue_df, n_snps > 1)
ggplot(data = tissue_df_nsnp_2, mapping = aes(x= reorder(tissue, -samples_eur))) +
  geom_bar() + 
  theme_bw() +
  guides(x = guide_axis(angle = 90)) +
  geom_text(aes(label = ..count..), stat = "count", vjust = -1, colour = "black", size = 2.5) +
  labs(x = "GTEx Tissue", y = paste("Number of ",QTL,"s \n across all genes",sep=""))
if (source == "predixcan") {
  ggsave(paste("./figures/summary_results/total_",QTL,"s_across_tissues_with_nsnps_more_than_1_predixcan_results.png",sep=""), dpi = 300)
} else if (source == "dbtables") {
  ggsave(paste("./figures/summary_results/total_",QTL,"s_across_tissues_with_nsnps_more_than_1_dbtables.png",sep=""), dpi = 300)
}
