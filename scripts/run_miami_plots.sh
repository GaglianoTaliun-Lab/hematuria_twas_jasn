#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=0:30:00
#SBATCH --job-name=miami_plots
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=2G

module load StdEnv/2020
module load r/4.1.0

export R_LIBS=~/.local/R/$EBVERSIONR/.

# need to indicate at least 2 arguments:
# 1) sQTL or eQTL
# 2) tissue index (1-49) # Note: 30 is for Kidney Cortex, 40 is for Skin(not sun exp), 34 is for skeletal muscle  and 49 for Whole Blood
# 3-X) one argument per gene to annotate in the TWAS

Rscript create_miami_plot.R sQTL 34
