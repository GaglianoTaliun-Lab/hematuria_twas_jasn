# Example to generate a Manhattan plot using the ggplot2 package.
# The function highlights SNPs that surpass the indicated threshold and draws a horizontal line on the threshold.

# Import libraries ---------------------------------------------------------------
library(here)
library(dplyr)
library(ggplot2)
library(readr)
library(readr)
library(ggrepel)
library(RColorBrewer)
library(stringr)

# Arguments -----------------------------------------------------------------------
project_dir = "/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project"
data_dir = stringr::str_c(project_dir, "/spredixcan_output/data_miami_plots")
multi_dir = stringr::str_c(project_dir,"/multixcan_output/data_miami_plots")

args <- commandArgs(TRUE)
QTL_type <- as.character(args[1]) # GWAS or sQTL or eQTL
twas_type <- as.character(args[2]) # GWAS or multitissue or kidney

# Import the function:
source(here(project_dir, "scripts", "ggmanhattan_function.R"))

# Suffix name for the output manhattan plot:
if (QTL_type == "GWAS") {
  name_plot = QTL_type
} else { name_plot = paste(QTL_type,twas_type,sep="_") }

out_dir = stringr::str_c(project_dir, "/spredixcan_output/figures")
out_png = stringr::str_c("manhattan_", name_plot, ".png")
out_pdf = stringr::str_c("manhattan_", name_plot, ".pdf")

# Main ----------------------------------------------------------------------------

# read TWAS/GWAS results:
if (QTL_type == "eQTL" & twas_type == "kidney") {
  plot_file <- read.table(here(data_dir, "hematuria_TOPMed_imputed_eqtl_Kidney_Cortex.tsv"), header = T, sep = "\t", colClasses = c("character","integer", "numeric", "numeric")) %>%
  filter(., chromosome.x %in% c(1:22)) %>%
  filter(., !is.na(start_position))
  data_color = "purple"
} else if (QTL_type == "sQTL" & twas_type == "kidney") {
  plot_file <- read.table(here(data_dir, "hematuria_TOPMed_imputed_sqtl_Kidney_Cortex.tsv"), header = T, sep = "\t", colClasses = c("character","integer", "numeric", "numeric"))
  data_color = "orange"
} else if (QTL_type == "eQTL" & twas_type == "multitissue") {
  plot_file <- read.table(here(multi_dir, "hematuria_eqtl_smultixcan_mapped.tsv"), header = T, sep = "\t", colClasses = c("character","integer", "numeric", "numeric")) %>%
  filter(., chromosome %in% c(1:22)) %>%
  filter(., !is.na(start_position))
  data_color = "purple"
} else if (QTL_type == "sQTL" & twas_type == "multitissue") {
  plot_file <- read.table(here(multi_dir, "hematuria_sqtl_smultixcan_mapped.tsv"), header = T, sep = "\t", colClasses = c("character","integer", "numeric", "numeric"))
  data_color = "orange"
} else if (QTL_type == "GWAS") {
  plot_file <- read.table(here(data_dir, "hematuria_GWAS_p_cutoff_0.1.tsv"), header = T, sep = "\t") %>%
    filter(., chrom != "X") %>%
    mutate(chr_pos = str_c(chrom, ":", pos)) %>%
    dplyr::select(chr_pos, chrom, pos, pval)
    dups <- plot_file[duplicated(plot_file),] %>% dplyr::select(chr_pos)
  plot_file <- plot_file %>% filter(., !chr_pos %in% dups)
  plot_file$chrom = as.numeric(plot_file$chrom)
  data_color = "red"
} else { stop("Need to choose either GWAS or eQTL or sQTL for the first argument.", call.=FALSE) }

# pvalue line thresholds:
if (QTL_type == "eQTL" & twas_type == "kidney") {
  pval_thres = 1.54e-5
} else if (QTL_type == "sQTL" & twas_type == "kidney") {
  pval_thres = 1.36e-5
} else if (QTL_type == "eQTL" & twas_type == "multitissue") {
  pval_thres = 3.73e-6
} else if (QTL_type == "sQTL" & twas_type == "multitissue") {
  pval_thres = 2.84e-6
} else if (QTL_type == "GWAS" & twas_type == "GWAS") {
  pval_thres = 5e-08
} else { stop("ERROR: Check arguments.", call.=FALSE) } 

# ylim thresholds:
if (QTL_type == "eQTL" & twas_type == "kidney") {
  ylim_thres = 20
} else if (QTL_type == "sQTL" & twas_type == "kidney") {
  ylim_thres = 25
} else if (QTL_type == "eQTL" & twas_type == "multitissue") {
  ylim_thres = 30
} else if (QTL_type == "sQTL" & twas_type == "multitissue") {
  ylim_thres = 25
} else if (QTL_type == "GWAS" & twas_type == "GWAS") {
  pval_thres = 50
} else { stop("ERROR: Check arguments.", call.=FALSE) }

# run function:
gg.manhattan(
  df_file = plot_file,
  threshold = pval_thres,
  ylims = c(0, ylim_thres),
  colours = c("#666666", "#333333"),
  highlight_color = data_color)

ggsave(here(out_dir,out_png), width = 50, height = 40, units = "cm")
ggsave(here(out_dir,out_pdf), width = 50, height = 40, units = "cm")
