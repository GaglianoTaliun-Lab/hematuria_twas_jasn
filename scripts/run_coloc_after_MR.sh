#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=1:00:00
#SBATCH --job-name=coloc
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G

module load StdEnv/2020
module load r/4.1.0

export R_LIBS=~/.local/R/$EBVERSIONR/.

# Rscript coloc_after_MR.R
Rscript coloc.R
