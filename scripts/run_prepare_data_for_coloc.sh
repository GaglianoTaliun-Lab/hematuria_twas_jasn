#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=6:00:00
#SBATCH --job-name=prepare_data_for_coloc
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G

module load StdEnv/2020
module load r/4.1.0

export R_LIBS=~/.local/R/$EBVERSIONR/.

Rscript prepare_data_for_coloc.R
