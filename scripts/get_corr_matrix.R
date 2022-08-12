# Description: get correlation matrix from LDSC results

# Load packages -----------------------------------------------------------

library(here)
library(tidyverse)
library(stringr)
library(data.table)
library(readr)

# Set arguments -----------------------------------------------------------

project_dir <- "/home/fridald4/projects/def-gsarah/fridald4/hematuria_project"

# phenotypes to include in the matrix
phenotypes_prox <- read.table(here(project_dir, "MR_analysis", "ldsc_corr", "list_traits.txt"), sep = "\t", header = F) %>% as.array() %>% sort()

# path for sumstats.gz:
gwas_dir <- here(project_dir, "MR_analysis", "ldsc_corr")

# path for rg output
out_dir <- here(project_dir,"MR_analysis", "ldsc_corr", "rg_output")

# number of lines to skip in the ldsc output:
n_skip=60
# number of lines to read from the output:
n_read=1

###### Extracting LDSC results ######

# Extract all outputs into single file
file_paths <-
  list.files(
    out_dir,
    pattern = "^rg_",
    full.names = T
  )

files <-
  vector(mode = "list",
         length = length(file_paths))

for(i in 1:length(file_paths)){
  
  files[[i]] <-
    read_table(
      file = file_paths[i],
      skip = n_skip, # Number of lines to skip in log file
      n_max = n_read # Number of lines to read
    )
  
}

# extract rg values in matrix
all_rg <-
  files %>%
  rbindlist(., fill = TRUE) %>%
  dplyr::mutate(
    p1 = basename(p1) %>%
      stringr::str_remove(".sumstats.gz") %>%
      stringr::str_remove(str_c(gwas_dir,"/")),
    p2 = basename(p2) %>%
      stringr::str_remove(".sumstats.gz") %>%
      stringr::str_remove(str_c(gwas_dir,"/"))
  ) %>%
  select(
    p1,
    p2,
    rg,
    se,
    z,
    p,
    h2_obs,
    h2_obs_se,
    h2_int,
    h2_int_se,
    gcov_int,
    gcov_int_se
  )

# Save data ---------------------------------------------------------------

write.table(
  all_rg,
  file = here(out_dir, stringr::str_c("CORR_MATRIX_LDSC.txt")),
  sep = "\t",
  row.names = F,
  quote = F
)
