library(purrr)
library(dplyr)
library(ggplot2)

# set home directory
home_dir = "/home/fridald4/projects/def-gsarah/fridald4/hematuria_project/spredixcan_output/heatmap_data/"

QTL="sqtl"

if (QTL == "eqtl") {
	# read filelist of significant genes with pattern for eqtl:
	filelist_sign <- list.files(pattern = "^hematuria_TOPMed_imputed_eqtl_significant", path = home_dir, full.names=TRUE)
	# read actual files and save in list
	tbl <- lapply(filelist_sign, function(x)read.table(x, header=T, sep= "\t"))
	# join list elements into a dataframe by column content
	tbl_sign <- reduce(tbl, full_join, by = c("gene_name", "tissue", "pvalue","zscore","chromosome","position"))
} else {
	# read filelist of significant genes with pattern for sqtl:
        filelist_sign <- list.files(pattern = "^hematuria_TOPMed_imputed_sqtl_significant", path = home_dir, full.names=TRUE)
        # read actual files and save in list
        tbl <- lapply(filelist_sign, function(x)read.table(x, header=T, sep= "\t"))
        # join list elements into a dataframe by column content
        tbl_sign <- reduce(tbl, full_join, by = c("gene_name", "tissue", "pvalue","zscore","chromosome","position"))
}

# for eQTLs, only plot selected genes in chr6:
if (QTL == "eqtl") {
 	tbl_sign_v2 <- filter(tbl_sign, chromosome != 6 | gene_name == "HLA-B" | "gene_name" == "HLA-C" | "gene_name" == "TGFB1" | "gene_name" == "CCDC97")
}

# plot heatmap with ggplot
if (QTL == "eqtl") {data_plot = tbl_sign_v2} else {data_plot = tbl_sign}

heatmap <- ggplot(data = data_plot, mapping = aes(x = tissue, y = reorder(gene_name, c(chromosome, position)), fill = zscore)) +
  geom_tile(colour = "black", lwd = 0.5, linetype = 1) +
  xlab(label = "tissue") +
  theme_classic() +
  scale_fill_gradient2(low = "blueviolet",
                       mid = "#FFFFCC",
                       high = "orangered") +
  guides(fill = guide_colourbar(title = "Z-score"), x = guide_axis(angle = 90))

# save plot
dir_plot = "/home/fridald4/projects/def-gsarah/fridald4/hematuria_project/spredixcan_output/figures/"
ggsave(plot=heatmap, filename=paste(dir_plot,"heatmap_",QTL,"_zscore_significant_all_tissues.pdf", sep=""), width = 25, height = 35, units = "cm")

# save long format table
dir_out = "/home/fridald4/projects/def-gsarah/fridald4/hematuria_project/spredixcan_output/heatmap_data/"
write.table(tbl_sign, paste(dir_out,"significant_genes_",QTL,"_for_heatmap.out",sep=""), sep="\t", quote = FALSE, row.names = FALSE)
