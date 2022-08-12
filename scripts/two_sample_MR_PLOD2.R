# Description: Mendelian Randomization analysis of hematuria-related traits using GWAS datasets
# (UK Biobank, CKDgen). Selected instruments based on eQTLs from GTEx (selected based on most
# significant tissue in GTEx web browser and using all independent predictors in that tissue for
# that gene in the mashr models for TWAS in MetaXcan)

# # tutorial: https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load libraries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(TwoSampleMR)
library(MRPRESSO)
library(here)
library(tidyverse)
library(stringr)
library(Cairo)
library(gridExtra)
library(grid)
library(gtable)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# current date - for output tables name
date=str_remove_all(Sys.Date(), "-")

project_dir <- here("/home/fridald4/projects/def-gsarah/fridald4/hematuria_project/MR_analysis")

# check available outcomes from the public API to obtain the ID (do it once)
# ao <- available_outcomes()

# exposure SNPs
snps_original <- c("rs7641322","rs57655223","rs73148997")
# palindromic: rs57655223 (proxy: rs1170389)
# not in CKDGen: rs7641322 (proxy: rs1449444) and rs73148997 (proxy: rs4681295)

snps_ckdgen <- c("rs1449444","rs1170389","rs4681295")
snps_ukbb <- c("rs7641322","rs1170389","rs73148997")
snps_meta <- c("rs1449444","rs1170389","rs4681295")
snps_stanzick <-c("rs7641322","rs1170389","rs73148997")

# function for pretty table
find_cells <- function(table, row, col, name="core-fg"){
  l <- table$layout
  unlist(Map(function(r, c) which(((l$t-1) == r) & ((l$l-1) == c) & (l$name == name)), row, col))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load files (exposures and outcomes)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1) Read exposure file for PLOD2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (obtained through Graham:"~/hematuria_project/MR_analysis/"
# from the eQTLs for Tibial Nerve downloaded from GTEx (EUR) GCP data:

exposure_dat <- read_exposure_data(
  filename = here(project_dir,"MR_exposures",'Skeletal_Muscle_EUR_PLOD2_exposure.tsv'),
  sep = '\t',
  snp_col = 'SNP',
  beta_col = 'beta',
  se_col = 'se',
  effect_allele_col = 'effect_allele',
  phenotype_col = 'Phenotype',
  units_col = 'units',
  other_allele_col = 'other_allele',
  eaf_col = 'eaf',
  samplesize_col = 'samplesize',
  ncase_col = 'ncase',
  ncontrol_col = 'ncontrol',
  gene_col = 'gene',
  pval_col = 'pval'
)

# filter exposures for CKDGEN
exposure_dat_ckdgen <- exposure_dat %>%
  filter(SNP %in% snps_ckdgen)

# filter exposures for UKBB data
exposure_dat_ukbb <- exposure_dat %>%
  filter(SNP %in% snps_ukbb)

# filter exposures for meta-analysis in ukbb+ckdgen
exposure_dat_meta <- exposure_dat %>%
  filter(SNP %in% snps_meta)

# filter exposures for meta-analysis from Stanzick et al. 2021
exposure_dat_stanzick <- exposure_dat %>%
  filter(SNP %in% snps_stanzick)

# 2) Read outcomes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# a. From the GWAS API for CKDgen data (eGFR): ieu-a-1105 = all inds; ieu-a-1104 = excluding diabetes inds.
outcome_eGFR_CKDG <- extract_outcome_data(exposure_dat_ckdgen$SNP,
                                          c('ieu-a-1105','ieu-a-1104'),
                                          proxies = 1,
                                          rsq = 0.8,
                                          align_alleles = 1,
                                          palindromes = 1,
                                          maf_threshold = 0.3) %>%
  mutate(outcome =
           case_when(
             id.outcome == "ieu-a-1104" ~ "eGFRcrea (CKDGen) - no DM",
             id.outcome == "ieu-a-1105" ~ "eGFRcrea (CKDGen) - ALL",
           ))

# b. From the GWAS API for CKDgen data (ACR): ieu-a-1107 = all inds; ieu-a-1101 = excluding diabetes inds.
outcome_ACR_CKDG <- extract_outcome_data(exposure_dat_ckdgen$SNP,
                                         c('ieu-a-1107','ieu-a-1101'),
                                         proxies = 1,
                                         rsq = 0.8,
                                         align_alleles = 1,
                                         palindromes = 1,
                                         maf_threshold = 0.3) %>%
  mutate(outcome =
           case_when(
             id.outcome == "ieu-a-1101" ~ "ACR (CKDGen) - no DM",
             id.outcome == "ieu-a-1107" ~ "ACR (CKDGen) - ALL",
           ))

# c. From the GWAS API for CKDgen data (eGFRcys): ieu-a-1106 = all inds
outcome_eGFRcys_CKDG <- extract_outcome_data(exposure_dat_ckdgen$SNP,
                                             c('ieu-a-1106'),
                                             proxies = 1,
                                             rsq = 0.8,
                                             align_alleles = 1,
                                             palindromes = 1,
                                             maf_threshold = 0.3) %>%
  mutate(outcome = "eGFRcys (CKDGen) - ALL")

# d. From the hematuria TopMed-imputed results (stored in Graham server):
outcome_hematuria_UKBB <- read_outcome_data(
  filename = here(project_dir,"MR_outcomes","hematuria_TOPMed_imputed_PLOD2_outcome.tsv"),
  phenotype_col = "phenotype",
  chr_col = "chr",
  pos_col = "pos",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "EAF",
  pval_col = "pval",
  samplesize_col = "N"
)

# e. From the ACR UK Biobank results (stored in Graham server: "ACR_sumstats"):
outcome_ACR_UKBB <- read_outcome_data(
  filename = here(project_dir,"MR_outcomes","ACR_UKBB_PLOD2_outcome.tsv"),
  phenotype_col = "phenotype",
  chr_col = "chr",
  pos_col = "pos",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "EAF",
  pval_col = "pval",
  samplesize_col = "N"
)

# f. From UK Biobank, composite phenotype (ACR + hematuria; stored in Graham server: "composite"):
outcome_composite_UKBB <- read_outcome_data(
  filename = here(project_dir,"MR_outcomes","composite_hematuria_ACR_UKBB_PLOD2_outcome.tsv"),
  phenotype_col = "phenotype",
  chr_col = "chr",
  pos_col = "pos",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "EAF",
  pval_col = "pval",
  samplesize_col = "N"
)

# g. ACR meta-analysis between UKBB and CKDGEN performed by me in metal
outcome_meta_ACR_all <- read_outcome_data(
  filename = here(project_dir,"MR_outcomes","output_metaanalysis_ACR_P3H2_PLOD1_PLOD2_all_1.tbl"),
  sep = "\t",
  snp_col = "MarkerName",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  pval_col = "P-value",
) %>%
  mutate(outcome = "meta-analysis_ACR")

# h. ACR meta-analysis between UKBB and CKDGEN performed by me in metal (excluding non-diabetes patients)
outcome_meta_ACR_noDM <- read_outcome_data(
  filename = here(project_dir,"MR_outcomes","output_metaanalysis_ACR_P3H2_PLOD1_PLOD2_noDM_1.tbl"),
  sep = "\t",
  snp_col = "MarkerName",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  pval_col = "P-value",
) %>%
  mutate(outcome = "meta-analysis_ACR_noDM")

# i. eGFRcreatinine meta-analysis from Stanzick et al. 2021 EUR (UKBB + CKDGEN)
outcome_meta_stanzick2021_eGFRcrea <- read_outcome_data(
  filename = here(project_dir,"MR_outcomes","metaanalysis_stanzick2021_EUR_eGFRcreatinine_PLOD2_outcome.tsv"),
  phenotype_col = "phenotype",
  chr_col = "chr",
  pos_col = "pos",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "EAF",
  pval_col = "pval",
  samplesize_col = "N"
)

# j. eGFRcystatin meta-analysis from Stanzick et al. 2021 EUR (UKBB + CKDGEN)
outcome_meta_stanzick2021_eGFRcys <- read_outcome_data(
  filename = here(project_dir,"MR_outcomes","metaanalysis_stanzick2021_EUR_eGFRcystatin_PLOD2_outcome.tsv"),
  phenotype_col = "phenotype",
  chr_col = "chr",
  pos_col = "pos",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "EAF",
  pval_col = "pval",
  samplesize_col = "N"
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main MR analysis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1) Harmonise instruments and outcomes:
# (action 3 means that all palyndromic SNPs are removed from analysis)
# a.
dat_eGFR_CKDG <- harmonise_data(exposure_dat_ckdgen, outcome_eGFR_CKDG, action = 3)
# b.
dat_ACR_CKDG <- harmonise_data(exposure_dat_ckdgen, outcome_ACR_CKDG, action = 3)
# c.
dat_eGFRcys_CKDG <- harmonise_data(exposure_dat_ckdgen, outcome_eGFRcys_CKDG, action = 3)
# d.
dat_hematuria_UKBB <- harmonise_data(exposure_dat_ukbb, outcome_hematuria_UKBB, action = 3)
# e.
dat_ACR_UKBB <- harmonise_data(exposure_dat_ukbb, outcome_ACR_UKBB, action = 3)
# f.
dat_composite_UKBB <- harmonise_data(exposure_dat_ukbb, outcome_composite_UKBB, action = 3)
# g.
dat_meta_ACR_all <- harmonise_data(exposure_dat_meta, outcome_meta_ACR_all, action = 3)
# h.
dat_meta_ACR_noDM <- harmonise_data(exposure_dat_meta, outcome_meta_ACR_noDM, action = 3)
# i.
dat_meta_stanzick2021_eGFRcrea <- harmonise_data(exposure_dat_stanzick, outcome_meta_stanzick2021_eGFRcrea, action = 3)
# j.
dat_meta_stanzick2021_eGFRcys <- harmonise_data(exposure_dat_stanzick, outcome_meta_stanzick2021_eGFRcys, action = 3)

# 2) Perform the analysis and save results in one large table
# Egger regression and other methods need at least 3 instruments!

mr_results_all_outcomes <- rbind(mr(dat_eGFR_CKDG, method_list = c("mr_ivw","mr_egger_regression")),
                                 mr(dat_ACR_CKDG, method_list = c("mr_ivw","mr_egger_regression")),
                                 mr(dat_eGFRcys_CKDG, method_list = c("mr_ivw","mr_egger_regression")),
                                 mr(dat_hematuria_UKBB, method_list = c("mr_ivw","mr_egger_regression")),
                                 mr(dat_ACR_UKBB, method_list = c("mr_ivw","mr_egger_regression")),
                                 mr(dat_composite_UKBB, method_list = c("mr_ivw","mr_egger_regression")),
                                 mr(dat_meta_ACR_all, method_list = c("mr_ivw","mr_egger_regression")),
                                 # mr(dat_meta_ACR_noDM, method_list = c("mr_ivw","mr_egger_regression")),
                                 mr(dat_meta_stanzick2021_eGFRcrea, method_list = c("mr_ivw","mr_egger_regression")),
                                 mr(dat_meta_stanzick2021_eGFRcys, method_list = c("mr_ivw","mr_egger_regression"))
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Post-MR: heterogeneity and outlier tests
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1) Perform heterogeneity test for eGFR in the Stanzick et al. 2021 meta-analysis
mr_heterogeneity(dat_eGFRcys_CKDG)

# 2) Test for horizontal pleiotropy using MR-Egger intercept test
pleiotropy_test <- mr_pleiotropy_test(dat_eGFRcys_CKDG)

# 3) Check for outliers in the hematuria_UKBB outcome using MRPRESSO
# paper: https://www.nature.com/articles/s41588-018-0099-7
# how to: https://github.com/rondolab/MR-PRESSO
# MR-PRESSO has three components, including: 1) detection of pleiotropy (MR-PRESSO global test); 
#   2) correction of pleiotropy via outlier removal (MR-PRESSO outlier test); and 
#   3) testing of significant distortion in the causal estimate before and after MR-PRESSO correction 
#   (MR-PRESSO distortion test).
#   Note: if the global test is not significant, then the other two components will not be performed!
#   Note2: needs at least 3 instruments.

mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
          SdOutcome = "se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_eGFRcys_CKDG, 
          NbDistribution = 1000,  SignifThreshold = 0.05)

# 4) Leave one out analysis
LOO <- mr_leaveoneout(dat_eGFRcys_CKDG)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate scatter plots 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Scatter plot (lines with methods used in results)
pdf(here(project_dir,"figures","PLOD2_scatter_plots.pdf"))
mr_scatter_plot(mr_results_all_outcomes, dat_eGFR_CKDG)
mr_scatter_plot(mr_results_all_outcomes, dat_ACR_CKDG)
mr_scatter_plot(mr_results_all_outcomes, dat_eGFRcys_CKDG)
mr_scatter_plot(mr_results_all_outcomes, dat_hematuria_UKBB)
mr_scatter_plot(mr_results_all_outcomes, dat_ACR_UKBB)
mr_scatter_plot(mr_results_all_outcomes, dat_composite_UKBB)
mr_scatter_plot(mr_results_all_outcomes, dat_meta_ACR_all)
# mr_scatter_plot(mr_results_all_outcomes, dat_meta_ACR_noDM)
mr_scatter_plot(mr_results_all_outcomes, dat_meta_stanzick2021_eGFRcrea)
mr_scatter_plot(mr_results_all_outcomes, dat_meta_stanzick2021_eGFRcys)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save main results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.table(mr_results_all_outcomes, here(project_dir,"results",str_c("MR_results_all_outcomes_PLOD2_",date,".tsv")), sep = "\t", row.names = F, quote = F)
write.table(pleiotropy_test, here(project_dir,"results",str_c("MR_pleiotropy_test_all_P3H2_",date,".tsv")), sep = "\t", row.names = F, quote = F)
write.table(LOO, here(project_dir,"results",str_c("MR_leave_one_out_analysis_eGFRcys_CKDGen_PLOD2_",date,".tsv")), sep = "\t", row.names = F, quote = F)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate html files for traits of interest
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mr_report(dat = dat_eGFRcys_CKDG,
         output_path = here(project_dir,"results"),
         output_type = "html",
         author = "FLD",
         study = "eGFRcys (CKDGen) - PLOD2 gene")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Output table in nice format of all results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mr_results_all_outcomes <- mr_results_all_outcomes %>%
  mutate(data_source =
           case_when(
             str_extract(outcome, "CKDGen") == "CKDGen" ~ "CKDGen",
             str_extract(outcome, "UKBB") == "UKBB" ~ "UKBiobank",
             str_extract(outcome, "meta") == "meta" ~ "meta-analysis (UKBiobank + CKDGen)"
           ),
         phenotype =
           case_when(
             outcome == "eGFRcrea (CKDGen) - no DM" ~ "eGFRcrea (no DM)",
             outcome == "eGFRcrea (CKDGen) - ALL" ~ "eGFRcrea",
             str_extract(outcome, "eGFRcrea") == "eGFRcrea" ~ "eGFRcrea",
             str_extract(outcome, "eGFRcys") == "eGFRcys" ~ "eGFRcys",
             str_extract(outcome, "ACR") == "ACR" ~ "ACR",
             str_extract(outcome, "hematuria") == "hematuria" ~ "Hematuria phecode 593",
             str_extract(outcome, "composite") == "composite" ~ "Composite (hematuria + ACR)"
           )
  )

ordered.outcomes <- c("eGFRcrea","eGFRcrea (no DM)","eGFRcys","ACR","Hematuria phecode 593","Composite (hematuria + ACR)")

mr_results_all_outcomes <- mr_results_all_outcomes %>%
  mutate(phenotype =  factor(phenotype, levels = ordered.outcomes)) %>%
  arrange(phenotype)

# gtable_add_grob:
# t	= a numeric vector giving the top extent of the grobs
# l	= a numeric vector giving the left extent of the grobs
# b	= a numeric vector giving the bottom extent of the grobs
# r	= a numeric vector giving the right extent of the grobs

tt_min <- ttheme_minimal()
mr <- mr_results_all_outcomes %>%
  mutate(beta = round(b, digits = 3),
         se = round(se, digits = 3),
         pvalue = round(pval, digits = 3)) %>%
  select(
    phenotype,
    data_source,
    exposure,
    method,
    nsnp,
    beta,
    se,
    pvalue
  )

mr <- tableGrob(mr, theme = tt_min, rows = NULL)
# add frames to table:
mr <- gtable_add_grob(mr,
                      grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                      t = 2, b = nrow(mr), l = 1, r = ncol(mr))
mr <- gtable_add_grob(mr,
                      grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                      t = 1, b = nrow(mr), l = 1, r = ncol(mr))
# add dashed lines:
mr <- gtable_add_grob(mr,
                      grobs = segmentsGrob(
                        x0 = unit(0,"npc"),
                        y0 = unit(0,"npc"),
                        x1 = unit(1,"npc"),
                        y1 = unit(0,"npc"),
                        gp = gpar(lwd = 2.0, lty = 3)),
                      t = 7, b = 7, l = 1, r = 8)
mr <- gtable_add_grob(mr,
                      grobs = segmentsGrob(
                        x0 = unit(0,"npc"),
                        y0 = unit(0,"npc"),
                        x1 = unit(1,"npc"),
                        y1 = unit(0,"npc"),
                        gp = gpar(lwd = 2.0, lty = 3)),
                      t = 11, b = 11, l = 1, r = 8)
mr <- gtable_add_grob(mr,
                      grobs = segmentsGrob(
                        x0 = unit(0,"npc"),
                        y0 = unit(0,"npc"),
                        x1 = unit(1,"npc"),
                        y1 = unit(0,"npc"),
                        gp = gpar(lwd = 2.0, lty = 3)),
                      t = 19, b = 19, l = 1, r = 8)

# highlight individual rows:
modify_cells <- function(mr, ids, gp=gpar()){
  for(id in ids) mr$grobs[id][[1]][["gp"]] <- gp
  return(mr)
}

ids <- find_cells(mr, 7, c(0:9), "core-fg")
mr <- modify_cells(mr, ids, gpar(fontface="bold"))

# add title:
title <- textGrob("Mendelian Randomization - PLOD2 gene",gp=gpar(fontsize=18))
padding <- unit(5,"mm")

mr <- gtable_add_rows(
  mr, 
  heights = grobHeight(title) + padding,
  pos = 0)
mr_table <- gtable_add_grob(
  mr, 
  title, 
  1, 1, 1, ncol(mr))

grid.newpage()
pdf(here(project_dir, "results", str_c("PLOD2_MR_results_table_format_",date,".pdf")), height = 8, width = 14)
grid.draw(mr_table)
dev.off()


