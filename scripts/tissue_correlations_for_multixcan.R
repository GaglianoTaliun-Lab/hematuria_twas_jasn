library(dplyr)

setwd("/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project/multixcan")

expression <- read.table("sharing_expression.txt", header = T, sep = "\t")
rownames(expression) <- expression[,1]
expression <- expression[,-1]
kdnctx_exp <- data.frame(expression$KDNCTX)
rownames(kdnctx_exp) <- rownames(expression)

splicing <- read.table("sharing_splice_junctions.txt", header = T, sep = "\t")
rownames(splicing) <- splicing[,1]
splicing <- splicing[,-1]
kdnctx_spl <- data.frame(splicing$KDNCTX)
rownames(kdnctx_spl) <- rownames(splicing)

eqtl <- read.table("sharing_cis_eqtls.mashrz.txt", header = T, sep = "\t")
rownames(eqtl) <- eqtl[,1]
eqtl <- eqtl[,-1]
kdnctx_eqtl <- data.frame(eqtl$KDNCTX)
rownames(kdnctx_eqtl) <- rownames(eqtl)

sqtl <- read.table("sharing_cis_sqtls.mashrz.txt", header = T, sep = "\t")
rownames(sqtl) <- sqtl[,1]
sqtl <- sqtl[,-1]
kdnctx_sqtl <- data.frame(sqtl$KDNCTX)
rownames(kdnctx_sqtl) <- rownames(sqtl)
