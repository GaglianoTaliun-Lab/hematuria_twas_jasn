#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=0:30:00
#SBATCH --job-name=metal
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2G

module load StdEnv/2020
module load gcc/9.3.0
module load metal/2011-03-25

HOMEPATH="/home/fridald4/projects/def-gsarah/fridald4/kidney_genetics_project"

metal ${HOMEPATH}/MR_analysis/meta_analysis_CKDgen_UKBB/metal_parameters_ACR_all.txt
metal ${HOMEPATH}/MR_analysis/meta_analysis_CKDgen_UKBB/metal_parameters_ACR_noDM.txt
