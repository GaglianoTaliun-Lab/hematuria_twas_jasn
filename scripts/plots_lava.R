library(here)
library(stringr)
library(dplyr)
library(tidyr)
library(tidygraph)
library(ggraph)
library(rtracklayer)
library(gghighlight)
library(forcats)

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/hematuria_project/MR_analysis/lava"
ref_dir = "/home/fridald4/projects/def-gsarah/fridald4/hematuria_project/reference_data"

# arguments
args <- commandArgs(TRUE)

# GWAS phenotypes that were included in the LAVA run:
datasets <- array(args[1:length(args)]) %>% sort()
lava_output_name <- str_c(datasets, collapse = ":")

bivar_all <- read.table(here(project_dir,"results","phenotypes_all_results_bivar.tsv"), sep = "\t", header = T)
bivar_bonf <- read.table(here(project_dir,"results","phenotypes_significant_results_bivar.tsv"), sep = "\t", header = T)

source(here(project_dir, "plots_RHR.R"))

# get number of bivariate tests
bivar_rds <- readRDS(here(project_dir,"results",stringr::str_c(lava_output_name,".bivar.lava.rds")))

# get number of bivariate tests
ntests = 0
for (i in 1:length(bivar_rds)) {
  condition = nrow(bivar_rds[[i]]) > 0
  if (condition == "TRUE" && length(condition) != 0) {
    ntests <- ntests + nrow(bivar_rds[[i]])
  }
}

# get pvalue bonferroni threshold
pvalue_bivar = 0.05/ntests

fct_disease = factor(datasets)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Heat map per LD block

bivar_all %>% 
  dplyr::filter(locus %in% unique(bivar_bonf$locus)) %>% 
  dplyr::mutate(
    rho_fill = 
      case_when(
        p < pvalue_bivar ~ round(rho, 2)
      ),
    locus_pos = str_c(chr,":",start,"-",stop)
  ) %>% 
  ggplot(
    aes(
      x = reorder(phen1, chr),
      y = reorder(phen2, chr),
      fill = rho_fill,
      label = round(rho, 2)
    )
  ) +
  geom_tile(colour = "black") +
  geom_text(
    size = 6
  ) +
  facet_wrap(vars(locus,locus_pos), ncol = 11) +
  scale_fill_distiller(palette = "RdBu", direction = -1, na.value = "#cccccc", limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) + theme_rhr +
  # theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        legend.key.width = unit(3, 'cm'),
        legend.key.height = unit(0.5, 'cm'),
        strip.text.x = element_text(size = 12, margin = margin(0.2,3,0.2,3, "cm"))) +
  labs(x = "eGFRcrea_stanzick2021", y = "eGFRcys_stanzick2021")

ggsave(here(project_dir, "figures", "phenotypes_heatmaps_per_LDblock.pdf"), width = 60, height = 55, units = "cm")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Annotation Plots per LD block

ref <- rtracklayer::import(here(ref_dir,"Homo_sapiens.GRCh37.87.gtf"))
ref <- ref %>% keepSeqlevels(c(1:22), pruning.mode = "coarse") 
ref <- ref[ref$type == "gene"]

loci_gr <-
  bivar_bonf %>%
  dplyr::count(locus, chr, start, stop, n_snps) %>%
  dplyr::arrange(locus) %>% 
  GenomicRanges::makeGRangesFromDataFrame(
    .,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    seqinfo = NULL,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "stop"
  )

fig_list = vector(mode = "list", length = length(loci_gr))

for(i in 1:length(loci_gr)){
  
  fig_list[[i]] <- 
    plot_locus(
    locus_gr = loci_gr[i], 
    ref = ref
    )
  
  names(fig_list)[i] <- str_c("locus_", loci_gr[i]$locus)
  
}

pdf(here(project_dir, "figures", "phenotypes_LDblock_annotations.pdf"))
fig_list
dev.off()

