#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=1:00:00
#SBATCH --job-name=get_nsnps_from_spredixcan
#SBATCH --output=slurm-%x-%a.out
#SBATCH --error=slurm-%x-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G

HOMEPATH="/home/fridald4/projects/def-gsarah/fridald4/hematuria_project/spredixcan_output/gene_mapped_results/"

awk '{OFS = "\t"} NR==1{print "tissue", $1, "n_snps"}' $HOMEPATH/hematuria_TOPMed_imputed_sqtl_Whole_Blood.tsv > /scratch/fridald4/snps_in_sites_per_tissue_from_spredixcan_results.tsv
awk '{OFS = "\t"} NR==1{print "tissue", $1, $2, "n_snps"}' $HOMEPATH/mapped_hematuria_TOPMed_imputed_eqtl_Whole_Blood.csv > /scratch/fridald4/snps_in_genes_per_tissue_from_spredixcan_results.tsv

# create new files to print only tissue, gene and n_snps_used:
for i in {1..49}
do
	TISSUE=$(sed -n ${i}p /scratch/fridald4/hematuria_scripts/list_of_tissues.txt)
	awk -v tissue="$TISSUE" '{OFS = "\t"} NR>1{print tissue, $1, $9}' $HOMEPATH/hematuria_TOPMed_imputed_sqtl_${TISSUE}.tsv >> /scratch/fridald4/snps_in_sites_per_tissue_from_spredixcan_results.tsv
	awk -v tissue="$TISSUE" '{OFS = "\t"} NR>1{print tissue, $1,$2, $7}' $HOMEPATH/mapped_hematuria_TOPMed_imputed_eqtl_${TISSUE}.csv >> /scratch/fridald4/snps_in_genes_per_tissue_from_spredixcan_results.tsv
done
