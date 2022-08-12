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
# 2-X) one argument per gene to annotate in the TWAS

Rscript create_miami_plot_multixcan_noFDR.R sQTL
