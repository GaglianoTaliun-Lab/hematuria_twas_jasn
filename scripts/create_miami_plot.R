# https://github.com/anastasia-lucas/hudson
# raw code: https://github.com/anastasia-lucas/hudson/blob/master/R/gmirror.R
# this code uses a modified version of the function named gmirror_v2.R 

# devtools::install_github('anastasia-lucas/hudson')

library(ggplot2)
library(gridExtra)
library(ggrepel)
library(hudson)

home_dir="/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project/spredixcan_output/data_miami_plots/"

# arguments
args <- commandArgs(TRUE)
QTL <- as.character(args[1])
i <- as.integer(args[2])

if (length(args) >= 3) {
	annotation_array <- array(args[3:length(args)])
} else {annotation_array = as.array(0)
}

# read TWAS results
if (QTL == "eQTL") {
  filelist <- list.files(pattern = "^hematuria_TOPMed_imputed_eqtl", path = home_dir, full.names=TRUE)
  twas <- lapply(filelist, function(x) read.table(x, header=T, sep= "\t"))
} else if (QTL == "sQTL") {
    filelist <- list.files(pattern = "^hematuria_TOPMed_imputed_sqtl", path = home_dir, full.names=TRUE)
    twas <- lapply(filelist, function(x) read.table(x, header=T, sep= "\t"))
  } else { stop("Need to choose either eQTL or sQTL for the first argument.", call.=FALSE) }

# read GWAS results
# gwas <- read.table(paste(home_dir, "hematuria_GWAS.tsv", sep=""), sep="\t", header = T)
gwas <- read.table(paste(home_dir, "hematuria_GWAS_p_cutoff_0.1.tsv", sep=""), sep="\t", header = T)
colnames(gwas) <- c("SNP","CHR","POS","pvalue")

# output directory
out_dir="/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project/spredixcan_output/figures/"

# read list of tissues
tissue <- read.table("/scratch/fridald4/hematuria_scripts/list_of_tissues.txt", header = F)

# read FDR thresholds
FDR <- read.csv("/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project/spredixcan_output/FDR_threshold_per_tissue.csv", header = T, sep = "\t")

# import modified function
source("/scratch/fridald4/hematuria_scripts/gmirror_v2.R")

# conditions to plot according to either eQTL or sQTL (pvalue cutoff for TWAS and FDR threshold):
if (QTL == "eQTL") {j = 2} else {j = 3}
if (QTL == "eQTL") {pval_twas = 1.60e-07} else {pval_twas = 1.38e-07}

twas_tmp <- as.data.frame(twas[[i]])
colnames(twas_tmp) <- c("SNP","CHR","POS","pvalue")
gmirror_v2(top=gwas, bottom=twas_tmp, tline=5e-08, bline=pval_twas, bline2= FDR[i,j],
	toptitle="Hematuria GWAS (ICD 593)", bottomtitle = paste("TWAS with GTEx - ",QTL,": ", tissue[i,1], sep=""), 
	highlight_p = c(5e-08,pval_twas), highlighter="purple", annotate_snp = annotation_array,
	file=paste(out_dir,"miami_plot_",QTL,"_",tissue[i,1],sep=""), type="png")

# to annotate SNPs/genes that passed a pvalue threshold, include in function: annotate_p = c(5e-08, 5.047e-08)
# to annotate particular genes/snps/sites, include in function: annotate_snps = c("rs123","gene1")
