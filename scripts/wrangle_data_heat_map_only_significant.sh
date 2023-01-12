#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=0:30:00
#SBATCH --array=1-49
#SBATCH --job-name=wrangle_for_heat_map
#SBATCH --output=slurm-%x-%a.out
#SBATCH --error=slurm-%x-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G

HOMEPATH="/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project"
TISSUE=$(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/fridald4/hematuria_scripts/list_of_tissues.txt)

# create new files to print only tissue, chr, pos, gene name, zscore and pvalues from annotated outputs:
awk -v tissue="$TISSUE" '{OFS = "\t"} NR==1{print "tissue", $2, "chromosome", "position", $3, $5} NR>1{print tissue, $2, $11, $12, $3, $5}' $HOMEPATH/spredixcan_output/gene_mapped_significant/hematuria_TOPMed_imputed_eqtl_${TISSUE}_significant.tsv > $HOMEPATH/spredixcan_output/heatmap_data/hematuria_TOPMed_imputed_eqtl_significant_${TISSUE}.out
awk -v tissue="$TISSUE" '{OFS = "\t"} NR==1{print "tissue", $2, "chromosome", "position", $3, $5} NR>1{print tissue, $1, $2, $3, $5, $7}' $HOMEPATH/spredixcan_output/sqtl_significant/hematuria_TOPMed_imputed_sqtl_${TISSUE}_significant.tsv > $HOMEPATH/spredixcan_output/heatmap_data/hematuria_TOPMed_imputed_sqtl_significant_${TISSUE}.out

# print only gene name and Z-score:
# awk -v tissue="$TISSUE" '{FS = OFS = "\t"} NR==1{print $2,tissue} NR>1{print $2, $3}' $HOMEPATH/spredixcan_output/eqtl/hematuria_TOPMed_imputed_eqtl_${TISSUE}.csv > $HOMEPATH/spredixcan_output/heatmap_data/hematuria_TOPMed_imputed_eqtl_${TISSUE}_zscore.out
# awk -v tissue="$TISSUE" '{FS = OFS = "\t"} NR==1{print $2,tissue} NR>1{print $2, $3}' $HOMEPATH/spredixcan_output/sqtl/hematuria_TOPMed_imputed_sqtl_${TISSUE}.csv > $HOMEPATH/spredixcan_output/heatmap_data/hematuria_TOPMed_imputed_sqtl_${TISSUE}_zscore.out

# print only gene name and pvalue:
# awk -v tissue="$TISSUE" '{FS = OFS = "\t"} NR==1{print $2,tissue} NR>1{print $2, $5}' $HOMEPATH/spredixcan_output/eqtl/hematuria_TOPMed_imputed_eqtl_${TISSUE}.csv > $HOMEPATH/spredixcan_output/heatmap_data/hematuria_TOPMed_imputed_eqtl_${TISSUE}_pval.out
# awk -v tissue="$TISSUE" '{FS = OFS = "\t"} NR==1{print $2,tissue} NR>1{print $2, $5}' $HOMEPATH/spredixcan_output/sqtl/hematuria_TOPMed_imputed_sqtl_${TISSUE}.csv > $HOMEPATH/spredixcan_output/heatmap_data/hematuria_TOPMed_imputed_sqtl_${TISSUE}_pval.out
