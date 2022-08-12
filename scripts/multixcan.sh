#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=2:00:00
#SBATCH --job-name=multixcan
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=5G

module load StdEnv/2020
module load python/3.7.9
module load scipy-stack/2021a

HOMEPATH="/home/fridald4/projects/def-gsarah/fridald4/renal_genetics_project"
QTL="sqtl" # or eqtl

python $HOMEPATH/MetaXcan/software/SMulTiXcan.py \
	--gwas_file $HOMEPATH/summary_statistics/hematuria_sumstats_with_chrpos_ID.txt.gz \
	--snp_column chrpos_ID \
	--models_folder $HOMEPATH/multixcan/MASHR/${QTL} \
	--models_name_pattern "mashr_(.*).db" \
	--snp_covariance $HOMEPATH/multixcan/gtex_v8_${QTL}_mashr_snp_covariance.txt.gz \
	--metaxcan_folder $HOMEPATH/multixcan/spredixcan_output/${QTL} \
	--metaxcan_filter "hematuria_TOPMed_imputed_(.*).csv" \
	--metaxcan_file_name_parse_pattern "(.*)_${QTL}_(.*).csv" \
	--effect_allele_column alt \
	--non_effect_allele_column ref \
	--beta_column beta \
	--se_column sebeta \
	--pvalue_column pval \
	--keep_non_rsid \
	--model_db_snp_key varID \
	--cutoff_condition_number 30 \
	--output $HOMEPATH/multixcan_output/hematuria_${QTL}_smultixcan.txt \
	--throw \
	--verbosity 7
