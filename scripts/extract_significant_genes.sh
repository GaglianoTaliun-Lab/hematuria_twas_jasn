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

HOMEPATH="/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project"
TISSUE=$(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/fridald4/hematuria_scripts/list_of_tissues.txt)

pval_genes=1.60e-07
pval_sites=1.38e-07

# sed -i 's/,/\t/g' $HOMEPATH/spredixcan_output/eqtl/hematuria_TOPMed_imputed_eqtl_${TISSUE}.csv
# sed -i 's/,/\t/g' $HOMEPATH/spredixcan_output/sqtl/hematuria_TOPMed_imputed_sqtl_${TISSUE}.csv

awk -F '\t' -v myvar=$pval_genes 'NR==1; NR>1 { if($5 <= myvar) print $1,$2,$3,$4,$5,$6,$10,$11,$12;}' $HOMEPATH/spredixcan_output/eqtl/hematuria_TOPMed_imputed_eqtl_${TISSUE}.csv > $HOMEPATH/spredixcan_output/eqtl_significant/hematuria_TOPMed_imputed_eqtl_${TISSUE}_significant.tsv
awk -F '\t' -v myvar=$pval_sites 'NR==1; NR>1 { if($5 <= myvar) print $1,$2,$3,$4,$5,$6,$10,$11,$12;}' $HOMEPATH/spredixcan_output/sqtl/hematuria_TOPMed_imputed_sqtl_${TISSUE}.csv > $HOMEPATH/spredixcan_output/sqtl_significant/hematuria_TOPMed_imputed_sqtl_${TISSUE}_significant.tsv

