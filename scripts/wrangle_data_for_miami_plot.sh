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

HOMEPATH="/home/fridald4/projects/def-gsarah/fridald4/renal_genetics_project"

TISSUE=$(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/fridald4/hematuria_scripts/list_of_tissues.txt)

# get data from TWAS:
awk '{FS = OFS = "\t"} {print $2,$11,$12,$5}' $HOMEPATH/spredixcan_output/gene_mapped_results/mapped_hematuria_TOPMed_imputed_eqtl_${TISSUE}.csv > $HOMEPATH/spredixcan_output/data_miami_plots/hematuria_TOPMed_imputed_eqtl_${TISSUE}.tsv

# get data from GWAS (this part of the script is not included in the array):
awk '{FS = OFS = "\t"} {print $5,$1,$2,$8}' $HOMEPATH/summary_statistics/hematuria_sumstats.tsv > $HOMEPATH/spredixcan_output/data_miami_plots/hematuria_GWAS.tsv
