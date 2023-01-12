library(dplyr)
library(tidyr)
library(arrow)

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project/GCP_all_associations/"

chr2_EUR_kidney <- read_parquet(paste0(project_dir,"Kidney_Cortex.v8.EUR.sqtl_allpairs.chr2.parquet"))
# chr2_all_kidney <- read_parquet(paste0(project_dir,"Kidney_Cortex.v8.cis_sqtl.all_pairs.chr2.parquet"))

# search for the COL4A4 exon 27 skipping event significant in kidney:
matches_EUR <- grep("chr2:227059623:227062530", chr2_EUR_kidney$phenotype_id)
chr2_EUR_kidney_exon27 <- chr2_EUR_kidney[matches,]
# matches_all <- grep("chr2:227059623:227062530", chr2_EUR_kidney$phenotype_id)
# chr2_all_kidney_exon27 <- chr2_all_kidney[matches,]

# search for COL4A4 SNP rs18898094 (chr2_227060058_C_T_b38) - to see if it is associated with other splicing events:
# matches_rsid_EUR <- grep("chr2_227060058_C_T_b38", chr2_EUR_kidney$variant_id)
# chr2_EUR_kidney_rs11898094 <- chr2_EUR_kidney[matches,]
# matches_rsid_all <- grep("chr2_227060058_C_T_b38", chr2_EUR_kidney$variant_id)
# chr2_all_kidney_rs11898094 <- chr2_all_kidney[matches,]

write.table(matches_EUR, paste0(project_dir,"matches_kidney_EUR_chr2:227059623:227062530.tsv"), sep = "\t", quote=F)
# write.table(matches_all, paste0(project_dir,"matches_kidney_all_chr2:227059623:227062530.tsv"), sep = "\t", quote=F)
# write.table(matches_rsid_EUR, paste0(project_dir,"matches_kidney_EUR_chr2_227060058_C_T_b3.tsv"), sep = "\t", quote=F)
# write.table(matches_rsid_all, paste0(project_dir,"matches_kidney_all_chr2_227060058_C_T_b3.tsv"), sep = "\t", quote=F)
