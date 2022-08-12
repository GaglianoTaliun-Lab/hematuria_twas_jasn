#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=0:30:00
#SBATCH --array=1-49
#SBATCH --job-name=extract
#SBATCH --output=slurm-%x-%a.out
#SBATCH --error=slurm-%x-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1G

HOMEPATH="/home/fridald4/projects/def-gsarah/fridald4/hematuria_project"
TISSUE=$(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/fridald4/hematuria_scripts/list_of_tissues.txt)

pval_genes=1.60e-07

awk -F '\t' -v myvar=$pval_genes 'NR==1; NR>1 { if($5 <= myvar) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13;}' $HOMEPATH/spredixcan_output/gene_mapped_results/mapped_hematuria_TOPMed_imputed_eqtl_${TISSUE}.csv > $HOMEPATH/spredixcan_output/gene_mapped_significant/hematuria_TOPMed_imputed_eqtl_${TISSUE}_significant.tsv

