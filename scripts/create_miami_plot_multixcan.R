# https://github.com/anastasia-lucas/hudson
# raw code: https://github.com/anastasia-lucas/hudson/blob/master/R/gmirror.R
# this code uses a modified version of the function named gmirror_v2.R 

# devtools::install_github('anastasia-lucas/hudson')

library(ggplot2)
library(gridExtra)
library(ggrepel)
library(hudson)

home_dir="/home/fridald4/projects/def-gsarah/fridald4/renal_genetics_project/multixcan_output/data_miami_plots/"

# arguments
args <- commandArgs(TRUE)
QTL <- as.character(args[1])

if (length(args) >= 2) {
	annotation_array <- array(args[2:length(args)])
} else {annotation_array = as.array(0)
}

# read TWAS results
if (QTL == "eQTL") {
  twas <- read.table(paste0(home_dir,"hematuria_eqtl_smultixcan_mapped.tsv"), header=T, sep= "\t")
} else if (QTL == "sQTL") {
  twas <- read.table(paste0(home_dir,"hematuria_sqtl_smultixcan_mapped.tsv"), header=T, sep= "\t")
  } else { stop("Need to choose either eQTL or sQTL for the first argument.", call.=FALSE) }

# read GWAS results
gwas <- read.table(paste0(home_dir,"hematuria_GWAS_p_cutoff_0.1.tsv"), sep="\t", header = T)
colnames(gwas) <- c("SNP","CHR","POS","pvalue")

# output directory
out_dir="/home/fridald4/projects/def-gsarah/fridald4/renal_genetics_project/multixcan_output/figures/"

# import modified function
source("/scratch/fridald4/hematuria_scripts/gmirror_v2.R")

# conditions to plot according to either eQTL or sQTL (pvalue cutoff for TWAS and FDR threshold):
if (QTL == "eQTL") {pval_twas = 3.735e-06} else {pval_twas = 2.836e-06}

# FDR pvalue
if (QTL == "eQTL") {fdr_twas = 7.566805e-05} else {fdr_twas = 3.733626e-05}

colnames(twas) <- c("SNP","CHR","POS","pvalue")
gmirror_v2(top=gwas, bottom=twas, tline=5e-08, bline=pval_twas, bline2=fdr_twas,
	toptitle="Hematuria GWAS (ICD 593)", bottomtitle = paste0("TWAS with GTEx - ",QTL,": MultiXcan (5 tissues)"), 
	highlight_p = c(5e-08, pval_twas), highlighter="purple", annotate_snp = annotation_array,
	file=paste0(out_dir,"miami_plot_",QTL,"_multixcan"), type="png")

# to annotate SNPs/genes that passed a pvalue threshold, include in function: annotate_p = c(5e-08, 5.047e-08)
# to annotate particular genes/snps/sites, include in function: annotate_snps = c("rs123","gene1")
