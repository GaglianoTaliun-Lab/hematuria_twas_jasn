# script to extract the SNPs used in the MASHR models for a specific gene/splice site across all 49 tissues in the GTEx.

library(RSQLite)
library(dplyr)
QTL="sqtl"

setwd(paste("/home/fridald4/projects/def-gsarah/fridald4/hematuria_project/MASHR/", QTL, "/mashr/", sep=""))

# connect to db table:

tissues <- read.table("/scratch/fridald4/hematuria_scripts/list_of_tissues.txt")

for (i in 1:49) {
  call_tissue <- tissues[i, 1]
  filename <- paste("mashr_", call_tissue,".db",sep="")  # corresponds to an eqtl
  sqlite.driver <- dbDriver("SQLite")
  db <- dbConnect(sqlite.driver,
                dbname = filename)
  # extract content:
  dbListTables(db)
  extra_info <- dbReadTable(db,"extra")
  weights_info <- dbReadTable(db, "weights")
  filter(weights_info, gene == "intron_16_58040698_58041617") %>% .[,2] -> tmp  # change the name of gene/site in this line (ENSgene name)
  if (length(tmp) > 0) tissues[i,1:length(tmp)+1] <- tmp
  else tissues[i,1] == "NA"
}

out_dir = "/home/fridald4/projects/def-gsarah/fridald4/hematuria_project/"
chr16_sQTLs <- tissues
write.table(chr16_sQTLs, paste(out_dir, "intron_16_58040698_58041617_SNPs_in_MASHR_models.tsv", sep=""), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
