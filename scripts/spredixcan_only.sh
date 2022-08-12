#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=0:30:00
#SBATCH --array=1-49
#SBATCH --job-name=spredixcan_sqtl
#SBATCH --output=slurm-%x-%a.out
#SBATCH --error=slurm-%x-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1G

module load StdEnv/2020
module load python/3.7.9
module load scipy-stack/2021a

HOMEPATH="/home/fridald4/projects/def-gsarah/fridald4/renal_genetics_project"
TISSUE=$(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/fridald4/hematuria_scripts/list_of_tissues.txt)
QTL="sqtl" # or eqtl

python $HOMEPATH/MetaXcan/software/SPrediXcan.py \
	--gwas_file $HOMEPATH/summary_statistics/hematuria_sumstats_with_chrpos_ID.txt.gz \
	--snp_column chrpos_ID \
	--effect_allele_column alt \
	--non_effect_allele_column ref \
	--beta_column beta \
	--se_column sebeta \
	--pvalue_column pval \
	--model_db_path $HOMEPATH/MASHR/${QTL}/mashr/mashr_${TISSUE}.db \
	--covariance $HOMEPATH/MASHR/${QTL}/mashr/mashr_${TISSUE}.txt.gz \
	--overwrite \
	--throw \
	--keep_non_rsid \
	--model_db_snp_key varID \
	--output_file $HOMEPATH/spredixcan_output/${QTL}/hematuria_TOPMed_imputed_${QTL}_${TISSUE}.csv
