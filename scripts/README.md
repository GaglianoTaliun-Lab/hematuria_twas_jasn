# hematuria_project

Project to perform TWAS on hematuria GWAS (UKBB with TOPMed imputations) across all GTEx tissues with the MASHR models, using S-PrediXcan (using Béluga cluster in Compute Canada).

- GWAS downloaded from: https://pheweb.org/UKB-TOPMed/pheno/593

- MASHR eQTL and sQTL models: https://predictdb.org/post/2021/07/21/gtex-v8-models-on-eqtl-and-sqtl/

- S-PrediXcan: https://github.com/hakyimlab/MetaXcan/wiki/S-PrediXcan-Command-Line-Tutorial

Short description of files in the repository:

- **list_of_tissues.txt**: list of 49 tissues in the GTEx to loop on.
- **wrangle_sumstats_one_liners.txt**: shell commands to adapt GWAS summary statistics according to S-PrediXcan input requirements.
- **spredixcan_only.sh**: run TWAS across all tissues
- **python_requirement.txt**: virtual environment for harmonization and imputation prior to performing TWAS (not used in final results)
- **virtual_env_steps.txt**: example on how to create a virtual environment on Compute Canada
- **extract_significant_genes.sh**: extract lines in TWAS results that passed a pvalue threshold
- **wrangle_data_heat_map_sh**: script to adapt TWAS results to plot on a heatmap
- **wrangle_data_heat_map_only_significant.sh**: same as above, but retaining only those genes/sites that were significant
- **wrangle_data_for_miami_plot.sh**: extracts columns from S-PrediXcan output (after mapping chr and pos) needed for miami plot and same for the GWAS sumstats
- **create_miami_plot.R**: function for obtaining a Miami plot of TWAS results
- **miami_plot_function.R**: general function for Miami plots, used in the script above
- **create_heatmap_plots.R**: function to obtain a heatmap from TWAS results
- **count_number_of_genes_in_mashr.R**: script to obtain proportions of genes in TWAS compared to total genes in models
- **filter_gene_map_one_liners.txt**: a few commands to download NCBI's gene information for genes, with chr and start/end positions (https://ftp.ncbi.nih.gov/gene/DATA/)
- **get_FDR.R**: script to obtain Storey's FDR using the R package "qvalue", for all tissues separately
- **map_genes.R**: include chromosome and start/end positions in the S-prediXcan outputs by using the data downloaded from NCBI with the 'filter_gene_map_one_liners.txt' commands
- **summary_plots_hematuria.R**: generates basic bar plots using as input the output table from 'count_number_of_genes_in_mashr.R', and computes correlations between the number of samples in the GTEx tissues and the genes in the models and in the TWAS
- **run_miami_plots.R**:
- **extract_significant_genes_mapped.sh**:
- **add_chrpos_sqtl.sh**:
