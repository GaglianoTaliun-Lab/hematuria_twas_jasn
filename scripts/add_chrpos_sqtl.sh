#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=1:00:00
#SBATCH --array=1-49
#SBATCH --job-name=add_chrpos_sqtl
#SBATCH --output=slurm-%x-%a.out
#SBATCH --error=slurm-%x-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G

HOMEPATH="/home/fridald4/projects/def-gsarah/fridald4/hematuria_project/spredixcan_output"
TISSUE=$(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/fridald4/hematuria_scripts/list_of_tissues.txt)

# output mapped data:
sed 's/_/\t/g'  $HOMEPATH/sqtl/hematuria_TOPMed_imputed_sqtl_${TISSUE}.csv | awk '{FS = OFS = "\t"} NR>1' | awk '{FS = OFS = "\t"} NR==1 {print "gene", "chromosome", "start_pos", "end_pos", "zscore", "effect_size", "pvalue", "var_g", "n_snps_used", "n_snps_in_cov", "n_snps_in_model"}  NR>1 {print $1"_"$2"_"$3"_"$4, $2, $3, $4, $9, $10, $11, $12, $16, $17, $18}' > $HOMEPATH/gene_mapped_results/hematuria_TOPMed_imputed_sqtl_${TISSUE}.tsv

# output mapped data columns for miami plots:
awk '{FS = OFS = "\t"} {print $1, $2, $3, $7}' $HOMEPATH/gene_mapped_results/hematuria_TOPMed_imputed_sqtl_${TISSUE}.tsv > $HOMEPATH/data_miami_plots/hematuria_TOPMed_imputed_sqtl_${TISSUE}.tsv

pval_sites=1.38e-07
# output significant results from mapped data:
awk -F '\t' -v var1=$pval_sites 'NR==1; NR>1 { if($7 <= var1) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11;}' $HOMEPATH/gene_mapped_results/hematuria_TOPMed_imputed_sqtl_${TISSUE}.tsv > $HOMEPATH/gene_mapped_significant/hematuria_TOPMed_imputed_sqtl_${TISSUE}_significant.tsv

