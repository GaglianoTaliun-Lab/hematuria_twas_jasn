# short script to check basic info in the MASHR models (i.e. number of SNP predictors, extract a specific SNP predictor)
# this script was not used to calculate the number of total tests, since it uses the information
# from dbtables and not from the TWAS results. To analyse the n_snps in TWAS output, use the
# script 'plot_nsnps.R'

library(RSQLite)
library(dplyr)
library(ggplot2)

# define if to look up in sQTL or eQTL tables:
QTL="eqtl"

setwd(paste("/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project/MASHR/", QTL, "/mashr/", sep=""))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# loop to extract info about number of SNP predictors per gene across all tissues
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# define if to look up in sQTL or eQTL tables:
QTL="eqtl"

# read list of tissues
tissues <- read.table("/scratch/fridald4/hematuria_scripts/list_of_tissues.txt")

# Create empty list for output:
tissue_snps_all <- list()

for (i in 1:49) {
  call_tissue <- tissues[i, 1]
  filename <- paste("mashr_", call_tissue,".db",sep="")
  sqlite.driver <- dbDriver("SQLite")
  db <- dbConnect(sqlite.driver,
                  dbname = filename)
  # extract content on nsnps:
  dbListTables(db)
  extra_info <- dbReadTable(db,"extra")
  weights <- dbReadTable(db,"weights")
  tissue_snps_all[[i]] <- data.frame("tissue" = rep(tissues[i,1], nrow(extra_info)), 
                                     "gene" = extra_info[,1], "gene_name" = extra_info[,2], "n_snps" = extra_info[,4])
}

# bind all elements in list ina  single dataframe (long format for plot)
tissue_snps_df_all <- bind_rows(tissue_snps_all)

### import file with sample size per tissue:
tissues_sample_size <- read.csv("/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project/GTEx_number_samples.csv", header = T)[,c(1,3)]
tissue_snps_df_all <- left_join(tissue_snps_df_all, tissues_sample_size, by = "tissue")

# export tables:
write.table(tissue_snps_df_all, "/scratch/fridald4/snps_in_sites_per_tissue_49tissues.tsv", sep = "\t", quote = F, row.names = F)

# NOTE: plots were generated in local computer by using the output table using the script 'plot_nsnps.R'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# find a specific SNP predictor from all db tables
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

QTL="eqtl"
tissues <- read.table("/scratch/fridald4/hematuria_scripts/list_of_tissues.txt")
tissue_snp <- list()

# lopp to extract info about number of SNP predictors per gene across all tissues:
for (i in 1:49) {
  call_tissue <- tissues[i, 1]
  filename <- paste("mashr_", call_tissue,".db",sep="")
  sqlite.driver <- dbDriver("SQLite")
  db <- dbConnect(sqlite.driver,
                  dbname = filename)
  # extract content on nsnps:
  dbListTables(db)
  dbReadTable(db,"weights") %>% filter(., rsid == "rs58261427") -> tissue_snp[[i]]
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# lopp to extract info about a specific gene across all tissues
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# define if to look up in sQTL or eQTL tables:
QTL="eqtl"

# which gene to query?
my_gene <- "PLOD2"

# read list of tissues
tissues <- read.table("/scratch/fridald4/hematuria_scripts/list_of_tissues.txt")

# create empty list for output:
tissue_gene <- list()

for (i in 1:49) {
  call_tissue <- tissues[i, 1]
  filename <- paste("mashr_", call_tissue,".db",sep="")
  sqlite.driver <- dbDriver("SQLite")
  db <- dbConnect(sqlite.driver,
                  dbname = filename)
  # extract content on nsnps:
  dbListTables(db)
  gene_ensID <- dbReadTable(db, "extra") %>% filter(., genename == my_gene) %>% .[1,1]
  extra_info <- dbReadTable(db,"extra") %>% filter(., gene == gene_ensID)
  weights <- dbReadTable(db,"weights") %>% filter(., gene == gene_ensID)
  tissue_gene[[i]] <- data.frame("tissue" = rep(tissues[i,1], nrow(extra_info)),
                                 "gene" = extra_info$gene, "gene_name" = extra_info$genename, 
                                 "N_predictors_by_tissue" = extra_info$n.snps.in.model, "rsid" = weights$rsid)
}

# bind all elements in list ina  single dataframe
tissue_gene <- bind_rows(tissue_gene)

write.table(tissue_gene, 
            paste0("/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project/spredixcan_output/predictors_for_MR/gene_",my_gene,".tsv"),
            sep = "\t", quote = F, row.names = F)
