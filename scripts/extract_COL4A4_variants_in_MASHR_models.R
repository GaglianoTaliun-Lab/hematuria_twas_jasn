library(RSQLite)
library(dplyr)
QTL="sqtl"

setwd(paste("/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project/MASHR/", QTL, "/mashr/", sep=""))

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
  filter(weights_info, gene == "ENSG00000081052.11") %>% .[,2] -> tmp
  if (length(tmp) > 0) tissues[i,1:length(tmp)+1] <- tmp
  else tissues[i,1] == "NA"
}

out_dir = "/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project/"
COL4A4_eQTLs <- tissues
write.table(COL4A4_eQTLs, paste(out_dir, "COL4A4_variants_in_MASHR_models.csv", sep=""), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
