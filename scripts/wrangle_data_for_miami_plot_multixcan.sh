#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=0:30:00
#SBATCH --job-name=wrangle_for_miami_plot
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G

HOMEPATH="/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project"

# get data from TWAS for eQTLs:
awk '{FS = OFS = "\t"} {print $2,$20,$22,$3}' $HOMEPATH/multixcan_output/hematuria_eqtl_smultixcan_mapped.txt > $HOMEPATH/multixcan_output/data_miami_plots/hematuria_eqtl_smultixcan_mapped.tsv

# get data from TWAS for sQTLs:
awk '{FS = OFS = "\t"} {print $1,$2,$3,$5}' $HOMEPATH/multixcan_output/hematuria_sqtl_smultixcan_mapped.txt > $HOMEPATH/multixcan_output/data_miami_plots/hematuria_sqtl_smultixcan_mapped.tsv

